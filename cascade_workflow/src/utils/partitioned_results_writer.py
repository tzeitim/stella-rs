"""
Partitioned Parquet Writer for Simulation Results

Efficiently writes simulation results to partitioned parquet files instead of thousands of individual files.
Uses Polars for optimal performance and automatic partitioning by key simulation parameters.
"""

import polars as pl
from pathlib import Path
import logging
from typing import Dict, Any, List, Optional
import fcntl
import time
from dataclasses import dataclass
from threading import Lock
import os

logger = logging.getLogger(__name__)


@dataclass
class PartitionConfig:
    """Configuration for parquet partitioning"""
    partition_columns: List[str]
    output_dir: Path
    buffer_size: int = 1000
    compression: str = "snappy"
    
    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)


class PartitionedResultsWriter:
    """
    Writes simulation results to partitioned parquet files for optimal storage and querying.
    
    Instead of thousands of individual files, creates a directory structure like:
    results/
    ├── cas9_tier_001/
    │   ├── solver_nj/
    │   │   └── data.parquet
    │   └── solver_maxcut/
    │       └── data.parquet
    ├── cas9_tier_002/
    │   ├── solver_nj/
    │   │   └── data.parquet
    │   └── solver_maxcut/
    │       └── data.parquet
    
    This enables:
    - Fast filtering: df.filter(pl.col("tier") == 2) loads only relevant partitions
    - Efficient storage: Better compression across similar records
    - Parallel processing: Different partitions can be processed independently
    """
    
    def __init__(self, partition_config: PartitionConfig):
        self.config = partition_config
        self.buffer: List[Dict[str, Any]] = []
        self.lock = Lock()
        
        logger.info(f"Initialized partitioned writer with partitions: {partition_config.partition_columns}")
        logger.info(f"Output directory: {partition_config.output_dir}")
        
    def add_result(self, result: Dict[str, Any]) -> None:
        """Add a single result to the buffer"""
        with self.lock:
            self.buffer.append(result.copy())
            
            if len(self.buffer) >= self.config.buffer_size:
                self._flush_buffer()
    
    def add_results_batch(self, results: List[Dict[str, Any]]) -> None:
        """Add multiple results at once"""
        with self.lock:
            self.buffer.extend([r.copy() for r in results])
            
            if len(self.buffer) >= self.config.buffer_size:
                self._flush_buffer()
    
    def _flush_buffer(self) -> None:
        """Flush buffered results to partitioned parquet files"""
        if not self.buffer:
            return

        try:
            # Convert to DataFrame
            df = pl.DataFrame(self.buffer)

            # Use custom partitioning with underscore naming instead of Hive format
            # Group by partition columns and write separately
            if self.config.partition_columns:
                partition_groups = df.group_by(self.config.partition_columns)

                for group_key, group_df in partition_groups:
                    # Create directory path with underscore format
                    partition_path = self.config.output_dir

                    # Handle both single and multiple partition columns
                    if isinstance(group_key, tuple):
                        for i, col in enumerate(self.config.partition_columns):
                            value = group_key[i]
                            # Zero-pad numeric values to 3 digits
                            if isinstance(value, (int, float)):
                                value = str(int(value)).zfill(3)
                            partition_path = partition_path / f"{col}_{value}"
                    else:
                        # Single partition column
                        value = group_key
                        # Zero-pad numeric values to 3 digits
                        if isinstance(value, (int, float)):
                            value = str(int(value)).zfill(3)
                        partition_path = partition_path / f"{self.config.partition_columns[0]}_{value}"

                    # Create directory if it doesn't exist
                    partition_path.mkdir(parents=True, exist_ok=True)

                    # Write group to separate file
                    output_file = partition_path / f"data_{int(time.time())}_{len(group_df)}.parquet"
                    group_df.write_parquet(
                        str(output_file),
                        compression=self.config.compression
                    )
            else:
                # No partitioning, write directly
                output_file = self.config.output_dir / f"data_{int(time.time())}_{len(df)}.parquet"
                df.write_parquet(
                    str(output_file),
                    compression=self.config.compression
                )

            logger.info(f"Flushed {len(self.buffer)} results to partitioned parquet")
            self.buffer.clear()
            
        except Exception as e:
            logger.error(f"Failed to flush buffer: {e}")
            # Don't clear buffer on error to avoid data loss
            raise
    
    def flush(self) -> None:
        """Manually flush any remaining buffered results"""
        with self.lock:
            self._flush_buffer()
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.flush()


class ConcurrentPartitionedWriter:
    """
    Thread-safe version that handles concurrent writes from multiple workers.
    Uses file locking to ensure data integrity when multiple processes write simultaneously.
    """
    
    def __init__(self, partition_config: PartitionConfig):
        self.config = partition_config
        self.local_buffer: List[Dict[str, Any]] = []
        self.lock_file = self.config.output_dir / ".writer.lock"
        
    def add_result(self, result: Dict[str, Any]) -> None:
        """Add result with concurrent safety"""
        self.local_buffer.append(result.copy())
        
        if len(self.local_buffer) >= self.config.buffer_size:
            self._flush_with_lock()
    
    def _flush_with_lock(self) -> None:
        """Flush with file system lock for concurrent access"""
        if not self.local_buffer:
            return
            
        # Use file-based locking for cross-process safety
        with open(self.lock_file, 'w') as lock_f:
            try:
                fcntl.flock(lock_f.fileno(), fcntl.LOCK_EX)
                
                df = pl.DataFrame(self.local_buffer)
                
                # Check if partitioned data already exists
                existing_files = list(self.config.output_dir.rglob("*.parquet"))
                
                if existing_files:
                    # Read existing data and append
                    existing_df = pl.scan_parquet(str(self.config.output_dir / "**/*.parquet")).collect()
                    combined_df = pl.concat([existing_df, df])
                    
                    # Remove old partitioned structure and rewrite
                    import shutil
                    for partition_col in self.config.partition_columns:
                        partition_dirs = list(self.config.output_dir.glob(f"{partition_col}=*"))
                        for pdir in partition_dirs:
                            if pdir.is_dir():
                                shutil.rmtree(pdir)
                    
                    combined_df.write_parquet(
                        str(self.config.output_dir),
                        compression=self.config.compression,
                        partition_by=self.config.partition_columns,
                        partition_chunk_size_bytes=128 * 1024 * 1024
                    )
                else:
                    # First write
                    df.write_parquet(
                        str(self.config.output_dir),
                        compression=self.config.compression,
                        partition_by=self.config.partition_columns,
                        partition_chunk_size_bytes=128 * 1024 * 1024
                    )
                
                logger.info(f"Flushed {len(self.local_buffer)} results to partitioned storage")
                self.local_buffer.clear()
                
            finally:
                fcntl.flock(lock_f.fileno(), fcntl.LOCK_UN)
    
    def flush(self) -> None:
        """Flush any remaining buffered results"""
        self._flush_with_lock()
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.flush()


def create_simulation_partition_config(output_dir: Path, buffer_size: int = 1000) -> PartitionConfig:
    """Create standard partition configuration for simulation results"""
    return PartitionConfig(
        partition_columns=["cas9_tier", "solver"],  # Primary partitioning
        output_dir=output_dir,
        buffer_size=buffer_size,
        compression="snappy"
    )


def read_partitioned_results(results_dir: Path, filters: Optional[Dict[str, Any]] = None) -> pl.DataFrame:
    """
    Read partitioned simulation results with optional filtering.
    
    Args:
        results_dir: Directory containing partitioned parquet files
        filters: Dict of column filters, e.g., {"cas9_tier": [1, 2], "solver": ["nj"]}
    
    Returns:
        Polars DataFrame with filtered results
        
    Example:
        # Read only tier 1 and 2 results with NJ solver
        df = read_partitioned_results(
            Path("results/"),
            filters={"cas9_tier": [1, 2], "solver": ["nj"]}
        )
    """
    
    # Construct file pattern for partitioned data
    pattern = str(results_dir / "**/*.parquet")
    
    # Create lazy scan
    lazy_df = pl.scan_parquet(pattern)
    
    # Apply filters if provided
    if filters:
        for col, values in filters.items():
            if isinstance(values, (list, tuple)):
                lazy_df = lazy_df.filter(pl.col(col).is_in(values))
            else:
                lazy_df = lazy_df.filter(pl.col(col) == values)
    
    # Collect and return
    return lazy_df.collect()


def get_partition_info(results_dir: Path) -> Dict[str, List[str]]:
    """
    Get information about available partitions.
    
    Returns:
        Dict mapping partition columns to their unique values
    """
    
    partition_info = {}
    
    # Look for partition directories (format: column=value)
    for path in results_dir.iterdir():
        if path.is_dir() and "=" in path.name:
            col, val = path.name.split("=", 1)
            if col not in partition_info:
                partition_info[col] = []
            partition_info[col].append(val)
    
    # Sort values for consistent output
    for col in partition_info:
        partition_info[col] = sorted(set(partition_info[col]))
    
    return partition_info


# Example usage and testing
if __name__ == "__main__":
    import tempfile
    import numpy as np
    
    # Test data
    test_results = []
    for tier in [1, 2, 3, 4]:
        for solver in ["nj", "maxcut", "greedy"]:
            for instance in range(10):
                test_results.append({
                    "gt_instance_id": instance,
                    "cas9_simulation_id": instance % 3,
                    "reconstruction_num": instance % 5,
                    "cas9_tier": tier,
                    "solver": solver,
                    "RF_distance": np.random.random(),
                    "triplets_distance": np.random.random(),
                    "cPHS": np.random.random(),
                    "computation_time_seconds": np.random.random() * 100
                })
    
    # Test partitioned writing
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        config = create_simulation_partition_config(temp_path / "partitioned_results", buffer_size=50)
        
        with PartitionedResultsWriter(config) as writer:
            for result in test_results:
                writer.add_result(result)
        
        # Test reading
        df = read_partitioned_results(
            temp_path / "partitioned_results",
            filters={"cas9_tier": [1, 2], "solver": ["nj"]}
        )
        
        print(f"Filtered results: {len(df)} rows")
        print(f"Unique tiers: {df['cas9_tier'].unique().sort().to_list()}")
        print(f"Unique solvers: {df['solver'].unique().to_list()}")
        
        # Show partition info
        info = get_partition_info(temp_path / "partitioned_results")
        print(f"Partition info: {info}")