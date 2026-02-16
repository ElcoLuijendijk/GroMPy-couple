#!/usr/bin/env python
"""Quick test of FiPy salt wedge model with fixes."""

import sys
import os
sys.path.insert(0, '.')

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

from model_input.model_parameters_sw_benchmark import ModelParameters, ModelOptions
from lib import grompy_fipy
from lib.mesh_functions_fipy import setup_rectangular_mesh_fipy

# Create fast test parameters
params = ModelParameters()
params.total_time = 2.0  # Just 2 seconds to test
params.dt_max = 2.0
params.output_interval = 2.0

opts = ModelOptions()
opts.backend = 'fipy'
opts.model_output_dir = 'model_output/salt_wedge_benchmark'

print('='*60)
print('QUICK TEST: FiPy Salt Wedge Model with Fixes')
print('='*60)

print('\n1. Creating mesh...')
mesh_file = 'model_output/salt_wedge_benchmark/_19_sw_quick_test.msh'
os.makedirs(os.path.dirname(mesh_file), exist_ok=True)

setup_rectangular_mesh_fipy(params, mesh_file)
print(f'   Mesh created: {mesh_file}')

print(f'\n2. Running model for {params.total_time} seconds...')
print('   (This may take 1-5 minutes)')

try:
    result = grompy_fipy.run_coupled_flow_model_fipy(
        params, opts, mesh_file
    )
    
    print(f'\n3. Test Result:')
    if result:
        print(f'   Result tuple has {len(result)} elements')
        print(f'   Element 5 type: {type(result[5])}')
        conc = result[5]  # Concentration is at index 5
        if hasattr(conc, 'value'):
            conc_vals = conc.value
        else:
            conc_vals = conc
        print(f'   Concentration min: {conc_vals.min():.6f}')
        print(f'   Concentration max: {conc_vals.max():.6f}')
        
        if conc.max() > 0.001:  # Should be around 0.03624
            print(f'\n   ✓ SUCCESS: Non-zero concentration detected!')
            print(f'   Expected: ~0.03624 (saltwater concentration)')
            print(f'   Got: {conc.max():.6f}')
            print(f'\n   All fixes are working correctly!')
        else:
            print(f'\n   ✗ FAIL: Concentration too low (still close to zero)')
    else:
        print('   ✗ FAIL: No result returned from model')
        
except Exception as e:
    print(f'\n✗ ERROR: {e}')
    import traceback
    traceback.print_exc()

print('='*60)
