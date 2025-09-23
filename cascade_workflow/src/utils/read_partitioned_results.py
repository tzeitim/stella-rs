#!/usr/bin/env python3
"""
Utility script for reading and analyzing partitioned simulation results
"""

import argparse
from pathlib import Path
import polars as pl
from partitioned_results_writer import read_partitioned_results, get_partition_info
import sys


def main():
    parser = argparse.ArgumentParser(description="Read and analyze partitioned simulation results")
    parser.add_argument('--results_dir', type=Path, required=True,
                       help='Directory containing partitioned results')
    parser.add_argument('--tier', type=int, nargs='+', 
                       help='Filter by Cas9 tiers (e.g., --tier 1 2)')
    parser.add_argument('--solver', type=str, nargs='+',
                       help='Filter by solvers (e.g., --solver nj maxcut)')
    parser.add_argument('--output', type=Path,
                       help='Optional output file for filtered results')
    parser.add_argument('--summary', action='store_true',
                       help='Show summary statistics')
    parser.add_argument('--info', action='store_true',
                       help='Show partition information only')
    
    args = parser.parse_args()
    
    if not args.results_dir.exists():
        print(f"Error: Results directory does not exist: {args.results_dir}")
        sys.exit(1)
    
    # Show partition info
    if args.info:
        partition_info = get_partition_info(args.results_dir)
        print("Available partitions:")
        for col, values in partition_info.items():
            print(f"  {col}: {values}")
        return
    
    # Build filters
    filters = {}
    if args.tier:
        filters['cas9_tier'] = args.tier
    if args.solver:
        filters['solver'] = args.solver
    
    # Read data
    print(f"Reading partitioned results from: {args.results_dir}")
    if filters:
        print(f"Applying filters: {filters}")
    
    try:
        df = read_partitioned_results(args.results_dir, filters=filters)
        print(f"Loaded {len(df)} results")
        
        if len(df) == 0:
            print("No results found matching the criteria")
            return
        
        # Show basic info
        print(f"Columns: {df.columns}")
        print(f"Unique tiers: {sorted(df['cas9_tier'].unique().to_list())}")
        print(f"Unique solvers: {sorted(df['solver'].unique().to_list())}")
        
        # Show summary statistics
        if args.summary:
            print("\n=== Summary Statistics ===")
            
            numeric_cols = ['RF_distance', 'triplets_distance', 'cPHS', 'computation_time_seconds']
            available_cols = [col for col in numeric_cols if col in df.columns]
            
            if available_cols:
                summary = df.select(available_cols).describe()
                print(summary)
                
                # Group by tier and solver
                print("\n=== Results by Tier and Solver ===")
                grouped = df.group_by(['cas9_tier', 'solver']).agg([
                    pl.col('RF_distance').mean().alias('avg_rf_distance'),
                    pl.col('triplets_distance').mean().alias('avg_triplets_distance'),
                    pl.col('cPHS').mean().alias('avg_cphs'),
                    pl.len().alias('count')
                ]).sort(['cas9_tier', 'solver'])
                print(grouped)
        
        # Save filtered results
        if args.output:
            df.write_parquet(args.output)
            print(f"Filtered results saved to: {args.output}")
        
        # Show first few rows
        print(f"\nFirst 5 rows:")
        print(df.head())
        
    except Exception as e:
        print(f"Error reading results: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()