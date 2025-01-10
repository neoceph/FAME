![Documentation Status](https://readthedocs.org/projects/fame-ud/badge/?version=latest)
[![PyPI version](https://img.shields.io/badge/TestPyPI-0.0.4-blue)](https://test.pypi.org/project/FAME-UD/)

## Development Instructions: [here](docs/developerReadme/development.md) &  Documentation: [here](https://fame-ud.readthedocs.io/en/latest/)

## Overview
FAME is a simulation tool designed to model the laser powder bed fusion (LPBF) additive manufacturing process using the Finite Volume Method (FVM). This tool helps in understanding the thermal and mechanical behavior of materials during the LPBF process.

## Features
- **Thermal Simulation**: Models the heat distribution and cooling rates.
- **Mechanical Simulation**: Analyzes stress and deformation.
- **Material Properties**: Supports various materials with customizable properties.
- **User-Friendly Interface**: Easy to set up and run simulations.

## Installation
To install FAME you need to have anaconda installed then clone the repository and install the required dependencies:
```bash
git clone https://github.com/neoceph/FAME.git
cd FAME
conda env create -f environment_linux.yaml
```

## Usage
To run a simulation, use the following command:
```bash
fame --input your_config_file.yaml
```
Replace `your_config_file.yaml` with your specific configuration file.

## Configuration
The configuration file should include parameters such as:
- Laser power
- Scan speed
- Layer thickness
- Material properties

Example configuration:
```yaml
simulation:
  domain:
    size: 
      x: [0, 0.5]
    divisions:
      x: [5]
    area: 10e-3
  
  material:
    name: "Aluminum"
    properties:  # Ensure all properties are nested under "properties"
      density:
        baseValue: 2700
        method: "constant"
        referenceTemperature: 298.15
      specificHeat:
        baseValue: 900
        method: "constant"
        referenceTemperature: 298.15
      thermalConductivity:
        baseValue: 1000
        method: "constant"
        referenceTemperature: 298.15

  boundaryConditions:
    parameters:
      temperature:
        variableType: "scalar"
        convectionCoefficient: 15
        emmissivity: 0.85
        ambientTemperature: 298
    x:
      0:  
        - type: "temperature"
          value: 100
      0.5:
        - type: "temperature"
          value: 500

  solver:
    method: "bicgstab"
    tolerance: 1e-8
    maxIterations: 1000
    preconditioner: "jacobi"

  timeControl:
    steadyState: true  # Indicates that this is a steady-state problem

  visualization:
    path: "./results"
    variableName: "temperature_cell"
```

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request.

## License
This project is licensed under the MIT License.

## Contact
For any questions or issues, please contact aamin1@udayton.edu.
