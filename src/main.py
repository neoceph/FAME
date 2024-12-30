import yaml

def loadConfig(file_path):
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)

def apply_boundary_conditions(boundary_conditions):
    for direction, conditions in boundary_conditions.items():
        for coord, bc in conditions.items():
            print(f"Applying boundary condition at {direction} = {coord}")
            print(f" - Type: {bc['type']}")
            if 'value' in bc:
                print(f" - Value: {bc['value']} K")
            if 'pressure' in bc:
                print(f" - Pressure: {bc['pressure']} Pa")

def main():
    config = loadConfig('simulation_config.yaml')
    boundary_conditions = config['simulation']['boundary_conditions']
    apply_boundary_conditions(boundary_conditions)

if __name__ == "__main__":
    main()
