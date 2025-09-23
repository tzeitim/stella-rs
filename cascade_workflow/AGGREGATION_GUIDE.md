# Large-Scale Results Aggregation Guide

## Problem
When dealing with tens of thousands of parquet result files, traditional batch concatenation becomes inefficient and memory-intensive.

## Solutions Implemented

### 1. Real-time Incremental Aggregation

**Automatic during reconstruction** - Results are added to a consolidated file as they're generated:

```python
# Already integrated in reconstruction_worker.py
# Creates: shared_dir/consolidated_results.parquet
```

This approach eliminates the need for post-processing by building the consolidated file incrementally.

### 2. Streaming Batch Aggregation

**For existing files** - Process thousands of files in manageable batches:

```bash
# Basic streaming aggregation
python cascade_workflow/src/utils/aggregate_results.py \
    --results_dir oishared_ny/results/ \
    --output consolidated_results.parquet \
    --method streaming \
    --batch_size 100
```

### 3. Database-Based Aggregation (DuckDB)

**For SQL-based analysis** - Use DuckDB for efficient aggregation and querying:

```bash
# DuckDB approach (requires: pip install duckdb)
python cascade_workflow/src/utils/aggregate_results.py \
    --results_dir oishared_ny/results/ \
    --output results.duckdb \
    --method duckdb
```

### 4. CSV Streaming (Most Memory Efficient)

**For maximum compatibility** - Stream to CSV format:

```bash
# CSV streaming approach
python cascade_workflow/src/utils/aggregate_results.py \
    --results_dir oishared_ny/results/ \
    --output consolidated_results.csv \
    --method csv
```

## Performance Comparison

| Method | Memory Usage | Processing Speed | File Size | SQL Query Support |
|--------|-------------|------------------|-----------|-------------------|
| Real-time Incremental | Low | Fast | Compact | No |
| Streaming Batches | Medium | Medium | Compact | Limited |
| DuckDB | Low | Very Fast | Compact | Full SQL |
| CSV Streaming | Very Low | Fast | Large | Limited |

## Recommended Approach

**For new runs**: Use the **real-time incremental aggregation** (already integrated) - it builds `shared_dir/consolidated_results.parquet` automatically.

**For existing files**: Use **DuckDB method** for best performance and query capabilities:

```bash
pip install duckdb

python cascade_workflow/src/utils/aggregate_results.py \
    --results_dir oishared_ny/results/ \
    --output analysis_results.duckdb \
    --method duckdb

# Then query with SQL:
python -c "
import duckdb
conn = duckdb.connect('analysis_results.duckdb')
print(conn.execute('SELECT solver, tier, AVG(cPHS) as avg_cphs FROM consolidated_results GROUP BY solver, tier').fetchall())
"
```

## Advanced Usage

### Custom Analysis Pipeline

```python
import polars as pl
from cascade_workflow.src.utils.streaming_aggregator import StreamingResultsAggregator

# Create custom aggregator
aggregator = StreamingResultsAggregator(
    consolidated_path="custom_analysis.parquet",
    batch_size=50  # Adjust based on memory
)

# Process with custom filtering
aggregator.process_directory(
    results_dir="oishared_ny/results/",
    file_pattern="*_vanilla_*.parquet"  # Only vanilla solver results
)
```

### Memory-Optimized Processing

For very large datasets (100k+ files):

```bash
# Use smaller batch sizes to reduce memory usage
python cascade_workflow/src/utils/aggregate_results.py \
    --results_dir oishared_ny/results/ \
    --output consolidated.parquet \
    --method streaming \
    --batch_size 25
```

## File Locking

The streaming aggregator includes file locking for concurrent access safety, allowing multiple workers to append results simultaneously without conflicts.