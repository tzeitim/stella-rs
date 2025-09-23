#!/usr/bin/env python3
"""
Standalone script for aggregating tens of thousands of parquet result files
"""
import argparse
import logging
from pathlib import Path
import polars as pl
from streaming_aggregator import StreamingResultsAggregator
import sys
import time

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def aggregate_existing_results(results_dir: Path, output_path: Path, batch_size: int = 100):
    """
    Aggregate all existing parquet files in a directory using streaming approach
    """
    logger.info(f"Starting aggregation of results from {results_dir}")
    start_time = time.time()
    
    aggregator = StreamingResultsAggregator(
        consolidated_path=output_path,
        batch_size=batch_size
    )
    
    # Check if this looks like partitioned results directory structure
    if any((results_dir / f"cas9_tier={i}").exists() for i in range(1, 10)):
        logger.info("Detected partitioned results structure, using '*.parquet' pattern")
        file_pattern = "*.parquet"
    else:
        logger.info("Using default '*_metrics.parquet' pattern")
        file_pattern = "*_metrics.parquet"

    aggregator.process_directory(
        results_dir=results_dir,
        file_pattern=file_pattern
    )
    
    elapsed = time.time() - start_time
    logger.info(f"Aggregation completed in {elapsed:.2f} seconds")
    
    # Verify final result
    if output_path.exists():
        final_df = pl.scan_parquet(str(output_path)).collect()
        logger.info(f"Final consolidated file contains {len(final_df)} records")
        logger.info(f"Columns: {final_df.columns}")
        return final_df
    else:
        logger.error("Aggregation failed - no output file created")
        return None

def database_approach_demo(results_dir: Path, db_path: Path):
    """
    Alternative: Use DuckDB for efficient aggregation
    """
    try:
        import duckdb
        
        logger.info("Using DuckDB approach for aggregation")
        
        # Create/connect to DuckDB database
        conn = duckdb.connect(str(db_path))
        
        # Create table if it doesn't exist
        parquet_files = list(results_dir.rglob("*_metrics.parquet"))
        if parquet_files:
            sample_df = pl.scan_parquet(str(parquet_files[0])).collect()
            
            # Create table with schema from first file
            conn.register('sample_df', sample_df)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS consolidated_results AS 
                SELECT * FROM sample_df WHERE FALSE
            """)
            
            # Insert data from all parquet files in batches
            batch_size = 50
            for i in range(0, len(parquet_files), batch_size):
                batch = parquet_files[i:i + batch_size]
                logger.info(f"Processing batch {i//batch_size + 1}/{(len(parquet_files)-1)//batch_size + 1}")
                
                for pf in batch:
                    try:
                        # Read and insert
                        df = pl.scan_parquet(str(pf)).collect()
                        conn.register('batch_df', df)
                        conn.execute("INSERT INTO consolidated_results SELECT * FROM batch_df")
                    except Exception as e:
                        logger.warning(f"Failed to process {pf}: {e}")
            
            # Export to parquet
            result = conn.execute("SELECT COUNT(*) FROM consolidated_results").fetchone()[0]
            logger.info(f"DuckDB aggregation completed: {result} records")
            
            # Optional: Export back to parquet
            conn.execute(f"COPY consolidated_results TO '{results_dir}/consolidated_duckdb.parquet' (FORMAT 'parquet')")
            
        conn.close()
        
    except ImportError:
        logger.error("DuckDB not available. Install with: pip install duckdb")

def streaming_csv_approach(results_dir: Path, output_path: Path):
    """
    Alternative: Stream to CSV format which is more append-friendly
    """
    logger.info("Using streaming CSV approach")
    
    parquet_files = list(results_dir.rglob("*_metrics.parquet"))
    if not parquet_files:
        logger.error("No parquet files found")
        return
    
    # Get schema from first file
    first_df = pl.scan_parquet(str(parquet_files[0])).collect()
    columns = first_df.columns
    
    # Write header
    with open(output_path, 'w') as f:
        f.write(','.join(columns) + '\n')
    
    # Stream append data
    batch_size = 50
    for i in range(0, len(parquet_files), batch_size):
        batch = parquet_files[i:i + batch_size]
        logger.info(f"Processing CSV batch {i//batch_size + 1}/{(len(parquet_files)-1)//batch_size + 1}")
        
        batch_dfs = []
        for pf in batch:
            try:
                df = pl.scan_parquet(str(pf)).collect()
                batch_dfs.append(df)
            except Exception as e:
                logger.warning(f"Failed to read {pf}: {e}")
        
        if batch_dfs:
            combined_batch = pl.concat(batch_dfs)
            # Append to CSV (mode='a' for append)
            combined_batch.write_csv(str(output_path), has_header=False, append=True)
    
    logger.info(f"CSV aggregation completed: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Aggregate thousands of parquet result files")
    parser.add_argument('--results_dir', type=Path, required=True,
                       help='Directory containing parquet result files')
    parser.add_argument('--output', type=Path, required=True,
                       help='Output path for consolidated results')
    parser.add_argument('--method', choices=['streaming', 'duckdb', 'csv'], 
                       default='streaming', help='Aggregation method to use')
    parser.add_argument('--batch_size', type=int, default=100,
                       help='Batch size for processing (default: 100)')
    
    args = parser.parse_args()
    
    if not args.results_dir.exists():
        logger.error(f"Results directory does not exist: {args.results_dir}")
        sys.exit(1)
    
    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    if args.method == 'streaming':
        result = aggregate_existing_results(
            results_dir=args.results_dir,
            output_path=args.output,
            batch_size=args.batch_size
        )
        
        if result is not None:
            logger.info("Aggregation successful!")
        else:
            sys.exit(1)
            
    elif args.method == 'duckdb':
        db_path = args.output.with_suffix('.duckdb')
        database_approach_demo(args.results_dir, db_path)
        
    elif args.method == 'csv':
        csv_path = args.output.with_suffix('.csv')
        streaming_csv_approach(args.results_dir, csv_path)

if __name__ == '__main__':
    main()