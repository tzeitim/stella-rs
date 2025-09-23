"""
Streaming aggregator for handling tens of thousands of parquet result files
"""
import polars as pl
from pathlib import Path
import logging
from typing import Optional, List, Dict, Any
import fcntl
import time

class StreamingResultsAggregator:
    """
    Efficiently aggregates results from many parquet files using streaming approach
    """
    
    def __init__(self, consolidated_path: Path, batch_size: int = 100):
        self.consolidated_path = Path(consolidated_path)
        self.batch_size = batch_size
        self.logger = logging.getLogger(__name__)
        
        # Ensure directory exists
        self.consolidated_path.parent.mkdir(parents=True, exist_ok=True)
        
    def append_results_batch(self, parquet_files: List[Path], use_lock: bool = True):
        """
        Append a batch of parquet files to consolidated results
        Uses file locking for concurrent worker safety
        """
        if not parquet_files:
            return
            
        try:
            # Read and concat the batch
            batch_dfs = []
            for pf in parquet_files:
                try:
                    df = pl.scan_parquet(str(pf)).collect()
                    batch_dfs.append(df)
                except Exception as e:
                    self.logger.warning(f"Failed to read {pf}: {e}")
                    continue
            
            if not batch_dfs:
                return
                
            batch_df = pl.concat(batch_dfs)
            
            # Append to consolidated file with locking
            if use_lock:
                self._append_with_lock(batch_df)
            else:
                self._append_direct(batch_df)
                
            self.logger.info(f"Appended {len(batch_dfs)} files to {self.consolidated_path}")
            
        except Exception as e:
            self.logger.error(f"Failed to append batch: {e}")
            raise
    
    def _append_with_lock(self, df: pl.DataFrame):
        """Append with file locking for concurrent access"""
        lock_file = self.consolidated_path.with_suffix('.lock')
        
        with open(lock_file, 'w') as lock_f:
            try:
                # Acquire exclusive lock
                fcntl.flock(lock_f.fileno(), fcntl.LOCK_EX)
                self._append_direct(df)
            finally:
                fcntl.flock(lock_f.fileno(), fcntl.LOCK_UN)
    
    def _append_direct(self, df: pl.DataFrame):
        """Direct append without locking"""
        if self.consolidated_path.exists():
            # Read existing, concat, and write back
            existing_df = pl.scan_parquet(str(self.consolidated_path)).collect()
            combined_df = pl.concat([existing_df, df])
        else:
            combined_df = df
            
        # Write back to parquet
        combined_df.write_parquet(str(self.consolidated_path))
    
    def process_directory(self, results_dir: Path, file_pattern: str = "*_metrics.parquet"):
        """
        Process all matching files in directory with streaming approach
        """
        parquet_files = list(Path(results_dir).rglob(file_pattern))
        self.logger.info(f"Found {len(parquet_files)} files to process")
        
        # Process in batches
        for i in range(0, len(parquet_files), self.batch_size):
            batch = parquet_files[i:i + self.batch_size]
            self.logger.info(f"Processing batch {i//self.batch_size + 1}/{(len(parquet_files)-1)//self.batch_size + 1}")
            
            self.append_results_batch(batch)
            
        self.logger.info(f"Completed processing {len(parquet_files)} files")


class IncrementalResultsWriter:
    """
    Alternative: Write results incrementally as they're generated
    """
    
    def __init__(self, output_path: Path):
        self.output_path = Path(output_path)
        self.buffer = []
        self.buffer_size = 1000  # Buffer up to N records before writing
        
    def add_result(self, result_dict: Dict[str, Any]):
        """Add a single result to buffer"""
        self.buffer.append(result_dict)
        
        if len(self.buffer) >= self.buffer_size:
            self.flush_buffer()
    
    def flush_buffer(self):
        """Write buffered results to file"""
        if not self.buffer:
            return
            
        df = pl.DataFrame(self.buffer)
        
        if self.output_path.exists():
            existing_df = pl.scan_parquet(str(self.output_path)).collect()
            combined_df = pl.concat([existing_df, df])
        else:
            combined_df = df
            
        combined_df.write_parquet(str(self.output_path))
        self.buffer.clear()
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.flush_buffer()