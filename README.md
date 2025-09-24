# stella-rs
High-performance phylogenetic analysis extensions for Cassiopeia lineage tracing, implemented in Rust with Python bindings.

It makes use of phylo-rs [S Vijendran · 2025](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06234-w). [phylo-rs crate](https://crates.io/crates/phylo)

## Cascade Workflow Diagram Generator

A configurable tool for generating visual workflow diagrams of the cascade job system based on YAML configuration files.

### Features

- **Automatic Configuration Parsing**: Reads YAML config files and extracts job parameters
- **Dynamic Layout**: Automatically adjusts diagram based on:
  - Number of GT instances
  - Number of Cas9 tiers
  - Number of enabled solvers
  - Total job counts
- **Color-coded Workflow Levels**: Different colors for each processing stage
- **Configuration Summary**: Displays key parameters and job counts
- **Batch Processing**: Generate diagrams for multiple configurations at once

### Usage

#### Single Diagram Generation
```bash
# Generate diagram from config file
python generate_workflow_diagram.py cascade_config.yaml

# Custom output path
python generate_workflow_diagram.py cascade_config_debug.yaml -o debug_workflow.png

# Show diagram after generation
python generate_workflow_diagram.py config.yaml --show
```

#### Batch Generation
```bash
# Generate diagrams for all cascade config files
./generate_all_workflow_diagrams.sh
```

### Example Configurations and Outputs

The generator automatically adapts to different scales:

- **Debug config**: 1 instance × 1 tier × 2 solvers = 4 total jobs
- **Main config**: 1000 instances × 4 tiers × 4 solvers = 21,000 total jobs  
- **Medium config**: 13 instances × 4 tiers × 7 solvers = 429 total jobs

### Workflow Structure

The generated diagrams show the 3-level job cascade:

1. **Level 1**: GT tree generation using `master_gt_worker.py`
2. **Level 2**: Cas9 recording across fidelity tiers using `cas9_recording_worker.py`
3. **Level 3**: Phylogenetic reconstruction with multiple solvers using `reconstruction_worker.py`

### Key Components

- **Entry Point**: `launch_cascade.py` - Orchestrates the entire workflow
- **Monitoring**: `job_monitor.py` - Continuous progress tracking
- **Configuration**: YAML files defining job parameters and solver settings
- **Output**: Consolidated parquet files with reconstruction metrics

### Files

- `generate_workflow_diagram.py` - Main diagram generator script
- `generate_all_workflow_diagrams.sh` - Batch processing script
- `cascade_config*.yaml` - Configuration files for different scenarios
- `workflow_cascade_config*.png` - Generated diagram outputs



