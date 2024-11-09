# FAME: FVM Based Laser Powder Bed Fusion Additive Manufacturing Process Simulation

## building the package

`python -m build`

## uploading the package to the test package host using "twine"
`twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose`

## Overview
FAME is a simulation tool designed to model the laser powder bed fusion (LPBF) additive manufacturing process using the Finite Volume Method (FVM). This tool helps in understanding the thermal and mechanical behavior of materials during the LPBF process.

## Features
- **Thermal Simulation**: Models the heat distribution and cooling rates.
- **Mechanical Simulation**: Analyzes stress and deformation.
- **Material Properties**: Supports various materials with customizable properties.
- **User-Friendly Interface**: Easy to set up and run simulations.

## Installation
To install FAME, clone the repository and install the required dependencies:
```bash
git clone https://github.com/yourusername/FAME.git
cd FAME
pip install -r requirements.txt
```

## Usage
To run a simulation, use the following command:
```bash
python simulate.py --config your_config_file.json
```
Replace `your_config_file.json` with your specific configuration file.

## Configuration
The configuration file should include parameters such as:
- Laser power
- Scan speed
- Layer thickness
- Material properties

Example configuration:
```json
{
    "laser_power": 200,
    "scan_speed": 1000,
    "layer_thickness": 0.03,
    "material": "Ti-6Al-4V"
}
```

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request.

## License
This project is licensed under the MIT License.

## Contact
For any questions or issues, please contact [your email].
