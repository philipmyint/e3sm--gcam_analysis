#!/usr/bin/env python3
"""
Quick test script to validate JSON configuration without running full analysis.
This helps debug configuration issues before submitting SLURM jobs.
"""

import json
import sys
import os
from pathlib import Path

def test_json_config(json_file):
    """Test JSON configuration for common issues."""
    
    print(f"Testing configuration file: {json_file}")
    print("=" * 60)
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        print(f"✓ JSON file loaded successfully")
        print(f"✓ Found {len(data)} configuration entries")
    except Exception as e:
        print(f"✗ Failed to load JSON: {e}")
        return False
    
    all_valid = True
    
    for i, entry in enumerate(data):
        print(f"\n--- Entry {i+1} ---")
        
        # Check required fields
        required_fields = ['variables', 'netcdf_substrings', 'simulation_path', 'output_file']
        for field in required_fields:
            if field not in entry:
                print(f"✗ Missing required field: {field}")
                all_valid = False
            else:
                print(f"✓ Has {field}")
        
        if 'variables' in entry and 'netcdf_substrings' in entry:
            # Check length consistency
            var_len = len(entry['variables'])
            netcdf_len = len(entry['netcdf_substrings'])
            
            print(f"Variables groups: {var_len}")
            print(f"NetCDF substring groups: {netcdf_len}")
            
            if var_len == netcdf_len:
                print("✓ Lengths match!")
            else:
                print(f"✗ Length mismatch: variables={var_len}, netcdf_substrings={netcdf_len}")
                all_valid = False
            
            # Show content structure
            print(f"Variables structure: {entry['variables']}")
            print(f"NetCDF structure: {entry['netcdf_substrings']}")
        
        # Check if simulation path exists
        if 'simulation_path' in entry:
            sim_path = Path(entry['simulation_path'])
            if sim_path.exists():
                print(f"✓ Simulation path exists: {sim_path}")
                # Count .nc files
                nc_files = list(sim_path.glob("*.nc"))
                print(f"  Found {len(nc_files)} .nc files")
            else:
                print(f"⚠ Simulation path not accessible: {sim_path}")
        
        # Check output directory
        if 'output_file' in entry:
            output_path = Path(entry['output_file']).parent
            print(f"Output directory: {output_path}")
            if not output_path.exists():
                print(f"⚠ Output directory doesn't exist, will be created")
    
    print("\n" + "=" * 60)
    if all_valid:
        print("✓ Configuration appears valid for SLURM submission!")
    else:
        print("✗ Configuration has issues - fix before submitting")
    
    return all_valid

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python test_config.py <json_file>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    success = test_json_config(json_file)
    sys.exit(0 if success else 1)