# Cascade Workflow - Organized Structure

This directory contains the organized structure for the LSF Cascading Job System for Cas9 lineage tracing analysis.

## Directory Structure

```
cascade_workflow/
├── src/                      # Source code
│   ├── core/                 # Core launcher and main scripts
│   │   └── launch_cascade.py # Main entry point
│   ├── workers/              # Worker scripts for distributed processing
│   │   ├── master_gt_worker.py        # Master ground truth generation
│   │   ├── cas9_recording_worker.py   # Cas9 recording simulation
│   │   ├── reconstruction_worker.py   # Tree reconstruction
│   │   └── job_monitor.py            # Job monitoring utility
│   └── utils/                # Utility modules
│       ├── config_loader.py  # Configuration management
│       └── solver_config.py  # Solver configuration
├── lsf_templates/            # LSF job submission templates
│   ├── master_gt_job.lsf               # Master job template
│   ├── cas9_recording_job_template.lsf # Cas9 recording job template
│   └── reconstruction_job_template.lsf # Reconstruction job template
├── configs/                  # Configuration files
│   ├── production/           # Production-ready configurations
│   │   ├── cascade_config.yaml        # Default configuration
│   │   ├── cascade_config_medium.yaml # Medium-scale configuration
│   │   └── cascade_config_pilot.yaml  # Pilot configuration
│   └── examples/             # Example and test configurations
│       ├── cascade_config_test.yaml        # Test configuration
│       ├── cascade_config_debug.yaml       # Debug configuration
│       ├── cascade_config_instance_001.yaml # Instance example 1
│       └── cascade_config_instance_002.yaml # Instance example 2
├── docs/                     # Documentation
│   ├── LSF_CASCADING_JOBS_USER_MANUAL.md  # User manual
│   ├── LSF_DISTRIBUTED_STRATEGY.md        # Architecture documentation
│   └── cascade_workflow_diagram.png       # Workflow diagram
└── scripts/                  # Utility scripts (placeholder for future scripts)
```

## Quick Start

### 1. Basic Launch
```bash
cd cascade_workflow
python src/core/launch_cascade.py --config configs/production/cascade_config.yaml
```

### 2. Test Launch (Local)
```bash
python src/core/launch_cascade.py --config configs/examples/cascade_config_test.yaml --local_test
```

### 3. Setup Only (No Job Submission)
```bash
python src/core/launch_cascade.py --config configs/production/cascade_config.yaml --setup_only
```

### 4. Monitor Jobs
```bash
python src/workers/job_monitor.py --shared_dir shared_medium --continuous
```

## Workflow Components

### Core Components
- **launch_cascade.py**: Main entry point that orchestrates the entire workflow
- **config_loader.py**: Handles YAML configuration parsing and validation
- **solver_config.py**: Manages solver-specific configurations

### Worker Scripts
- **master_gt_worker.py**: Generates ground truth phylogenetic trees
- **cas9_recording_worker.py**: Simulates Cas9 recording on cells
- **reconstruction_worker.py**: Reconstructs trees from Cas9 data
- **job_monitor.py**: Monitors job progress and status

### LSF Templates
Pre-configured job submission templates with placeholder substitution for:
- Resource allocation (CPU, memory, runtime)
- Job dependencies
- Output directories
- Environment setup

### Configuration Files

#### Production Configs
- **cascade_config.yaml**: Standard production configuration
- **cascade_config_medium.yaml**: Medium-scale runs (balanced resources)
- **cascade_config_pilot.yaml**: Small pilot runs for testing

#### Example Configs
- **cascade_config_test.yaml**: Minimal test configuration
- **cascade_config_debug.yaml**: Debug mode with verbose logging
- **cascade_config_instance_*.yaml**: Instance-specific examples

## Configuration Parameters

Key configuration sections:
- `simulation`: Tree generation parameters
- `recording`: Cas9 recording parameters
- `solvers`: Available reconstruction algorithms
- `execution`: Runtime parameters (instances, parallelism)
- `output`: Output directories and formats

## Job Flow

1. **Master GT Job** → Generates ground truth tree
2. **Cas9 Recording Jobs** (4x) → Simulates recording with different parameters
3. **Reconstruction Jobs** (7x per recording) → Reconstructs trees using different solvers
4. **Total**: 1 + 4 + (4 × 7) = 33 jobs per cascade

## Documentation

- [User Manual](docs/LSF_CASCADING_JOBS_USER_MANUAL.md) - Complete usage guide
- [Architecture](docs/LSF_DISTRIBUTED_STRATEGY.md) - System design and strategy

## Notes

- Original files remain untouched in the parent directory
- This organized structure facilitates easier maintenance and deployment
- All paths in job templates use placeholders that are substituted at runtime