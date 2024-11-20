![Documentation Status](https://readthedocs.org/projects/fame-ud/badge/?version=latest)
[![PyPI version](https://img.shields.io/badge/TestPyPI-0.1.0-blue)](https://test.pypi.org/project/FAME-UD/)

# Generating Sphinx documentation
From the root directory run the following commands [Powershell in windows or bash in linux]
- `./docs/make clean`
- `sphinx-apidoc -o ./docs/source/ ./src`
- `./docs/make html`
- `./docs/make latexpdf` to generate pdf. The appropriate latex compiler must be installed and available.

# Generating Docstring using VSCode extension

If you are using VS code, the [Python Docstring](https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring) extension can be used to auto-generate a docstring snippet once a function/class has been written. If you want the extension to generate docstrings in `Sphinx` format, you must set the `"autoDocstring.docstringFormat": "sphinx"` setting, under File > Preferences > Settings.

Note that it is best to write the docstrings once you have fully defined the function/class, as then the extension will generate the full dosctring. If you make any changes to the code once a docstring is generated, you will have to manually go and update the affected docstrings.

# FAME: FVM Based Laser Powder Bed Fusion Additive Manufacturing Process Simulation

## restructred text live preview on vscode

- need the extension `pip install esbonio` and 'esbonio'. After that make sure python path is manually setup if esbonio is having difficulty finding the python interpreter. You can do that by going to 'File->Preference->Settings' and finding Esbonio>Server:Python Path

## restructured text cheatsheet

Details of rST is found [here](https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html)

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
