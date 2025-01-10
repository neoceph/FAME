import yaml
import os
import sys
import argparse

# Adjust the import path to locate finiteVolumeMethod.py
from fame.FVM.finiteVolumeMethod import FVM

def loadInput(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Configuration file not found: {file_path}")
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run FAME FVM Simulation")
    parser.add_argument(
        '--input', 
        type=str, 
        default='config.yaml', 
        help="Path to the YAML input file"
    )
    
    args = parser.parse_args()
    input_path = args.input
    
    try:
        config = loadInput(input_path)
    except FileNotFoundError as e:
        print(e)
        sys.exit(1)
    
    # Instantiate and run the FVM simulation
    fvm_simulation = FVM(config)
    fvm_simulation.simulate()

if __name__ == "__main__":
    main()
