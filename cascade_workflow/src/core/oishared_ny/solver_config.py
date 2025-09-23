#!/usr/bin/env python3
"""
Centralized solver configuration for the LSF cascading job system.

This module defines all available solvers and their configurations in one place
to ensure consistency across all worker files.
"""

import cassiopeia as cass

# Available solvers configuration
SOLVERS_CONFIG = {
    'nj': {
        'name': 'Neighbor Joining',
        'class': lambda: cass.solver.NeighborJoiningSolver(add_root=True),
        'enabled': True
    },
    'maxcut': {
        'name': 'MaxCut',
        'class': lambda: cass.solver.MaxCutSolver(),
        'enabled': True
    },
    'greedy': {
        'name': 'Vanilla Greedy',
        'class': lambda: cass.solver.VanillaGreedySolver(),
        'enabled': True
    },
    'vanilla': {
        'name': 'Vanilla Greedy (alias)',
        'class': lambda: cass.solver.VanillaGreedySolver(),
        'enabled': True
    },
    'spectral': {
        'name': 'Spectral',
        'class': lambda: cass.solver.SpectralSolver(),
        'enabled': True
    },
    'smj': {
        'name': 'Shared Mutation Joining',
        'class': lambda: cass.solver.SharedMutationJoiningSolver(),
        'enabled': True
    },
    'dmj': {
        'name': 'Distance',
        'class': lambda: cass.solver.DistanceSolver(),
        'enabled': True
    },
    'ilp': {
        'name': 'ILP',
        'class': lambda: cass.solver.ILPSolver(),
        'enabled': False  # Disabled by default
    }
}

# Get list of enabled solvers
SOLVERS = [solver for solver, config in SOLVERS_CONFIG.items() if config['enabled']]

def get_solver_class(solver_name: str):
    """Returns the solver class based on the solver name"""
    if solver_name not in SOLVERS_CONFIG:
        raise ValueError(f"Unknown solver: {solver_name}. Available solvers: {list(SOLVERS_CONFIG.keys())}")
    
    if not SOLVERS_CONFIG[solver_name]['enabled']:
        raise ValueError(f"Solver '{solver_name}' is disabled. Enable it in SOLVERS_CONFIG to use.")
    
    return SOLVERS_CONFIG[solver_name]['class']()

def get_enabled_solvers():
    """Returns list of enabled solver names"""
    return SOLVERS

def get_all_solvers():
    """Returns list of all solver names (enabled and disabled)"""
    return list(SOLVERS_CONFIG.keys())

def enable_solver(solver_name: str):
    """Enable a solver"""
    if solver_name in SOLVERS_CONFIG:
        SOLVERS_CONFIG[solver_name]['enabled'] = True
        # Update SOLVERS list
        global SOLVERS
        SOLVERS = [solver for solver, config in SOLVERS_CONFIG.items() if config['enabled']]
    else:
        raise ValueError(f"Unknown solver: {solver_name}")

def disable_solver(solver_name: str):
    """Disable a solver"""
    if solver_name in SOLVERS_CONFIG:
        SOLVERS_CONFIG[solver_name]['enabled'] = False
        # Update SOLVERS list
        global SOLVERS
        SOLVERS = [solver for solver, config in SOLVERS_CONFIG.items() if config['enabled']]
    else:
        raise ValueError(f"Unknown solver: {solver_name}")

def validate_requested_solvers(requested_solvers: list) -> tuple[bool, list]:
    """
    Validate that all requested solvers are supported by the workflow.

    Args:
        requested_solvers: List of solver names from config

    Returns:
        tuple: (is_valid, unsupported_solvers)
               is_valid is True if all solvers are supported,
               unsupported_solvers is list of unsupported solver names
    """
    supported_solvers = get_enabled_solvers()
    unsupported = [solver for solver in requested_solvers if solver not in supported_solvers]
    return len(unsupported) == 0, unsupported

def get_supported_solvers_info() -> dict:
    """
    Get detailed information about all supported solvers.

    Returns:
        dict: Mapping of solver name to info dict with 'name' and 'enabled' status
    """
    return {
        solver: {
            'name': config['name'],
            'enabled': config['enabled']
        }
        for solver, config in SOLVERS_CONFIG.items()
        if config['enabled']
    }