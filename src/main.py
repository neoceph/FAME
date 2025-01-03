import yaml
import os
import sys
import argparse

# Adjust the import path to locate finiteVolumeMethod.py
from src.FVM.finiteVolumeMethod import FVM

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
        default='../configs/simulation_config.yaml', 
        help="Path to the YAML input file"
    )
    
    args = parser.parse_args()
    input_path = os.path.abspath(os.path.join(os.path.dirname(__file__), args.input))
    
    try:
        config = loadInput(input_path)
    except FileNotFoundError as e:
        print(e)
        sys.exit(1)
    
    # Instantiate and run the FVM simulation
    fvm_simulation = FVM(config)
    fvm_simulation.run()

if __name__ == "__main__":
    main()
