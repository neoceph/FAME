# Steps to setup environment for development

**1.** Make sure you have [anaconda](https://www.anaconda.com/download) installed and working. Clone the repository from github 
```bash
git clone https://github.com/neoceph/FAME
```
**2.** Create appropriate python environment by issuing the following command
```bash
conda env create -f environment_linux.yaml
```
If there are errors, try
```bash
conda env create -f environment_linux.yaml --force
```
If environment is created and needs updating the packages list after dependancy install (i.e. petsc4py) try
```bash
conda env update --name FAME --file environment_linux.yaml
```
If PETSc installation is failing make sure you have BLAS and LAPACK available

```bash
conda install -c conda-forge blas lapack
```
**3.** You are ready to contribute code under `.\fame` or `.\tests`.

### When its time to push changes read the instructions from [here](./pushingChanges.md)

## Download and installation of vtk from source if needed (Windows)

**1.** Make sure you have cmake installed and the path is added to environment variable.
```powershell
git clone https://gitlab.kitware.com/vtk/vtk.git
cd vtk
```
**2.** Create a directory titled build
```powershell
mkdir build
cd build
```
**3.** Configure the project using CMake. Make note of the appropriate 'Visual Studio' version specification. You also need to have write access to the `D:/destination/to/installation/of/VTK` directory. You might need to open command prompt or powershell as an administrator. Make sure you have <span style="color: red;">hdf5=YES.</span> set for installation setup.
```powershell
cmake .. -G "Visual Studio 17 2022" -A x64 `
  -DBUILD_SHARED_LIBS=ON `
  -DVTK_WRAP_PYTHON=ON `
  -DVTK_GROUP_ENABLE_Qt=NO `
  -DVTK_USE_SYSTEM_HDF5=OFF `
  -DVTK_MODULE_ENABLE_VTK_hdf5=YES `
  -DVTK_PYTHON_VERSION=3 `
  -DVTK_INSTALL_PYTHON_MODULE_DIR="$env:CONDA_PREFIX\Lib\site-packages" `
  -DPYTHON_EXECUTABLE="$env:CONDA_PREFIX\python.exe" `
  -DCMAKE_INSTALL_PREFIX="D:/destination/to/installation/of/VTK"
```
**4.** Install with multiple cpus. If you are using single cpu, it will take days.
```powershell
cmake --build . --config Release --target INSTALL --parallel
```
**5.** Once installation complete, create a `set_pythonpath.ps1` file for powershell or `set_pythonpath.bat` for CMD under directory `anaconda\installation\directory\etc\conda\activate.d`and set the VTK installation path.

If specific environment other than the `base` environment is of interest, the files need to be placed under `%CONDA_PREFIX%\etc\conda\deactivate.d` and `%CONDA_PREFIX%\etc\conda\activate.d`

powershell
```powershell
$env:PYTHONPATH = "D:/destination/to/installation/of/VTK/Lib/site-packages;" + $env:PYTHONPATH
```

CMD
```cmd
set PYTHONPATH=D:\ProgramFiles\VTK\Lib\site-packages;%PYTHONPATH%
```
**6.** You can check the PYTHONPATH variable with
powershell
```powershell
echo $env:PYTHONPATH
```
CMD
```CMD
echo %PYTHONPATH%
```

**7.** Setup clearing of PYTHONPATH for conda environment deactivation by creating `unset_pythonpath.ps1` or `unset_pythonpath.bat` with the following content

powershell
```powershell
$env:PYTHONPATH = ""
```

CMD
```cmd
@echo off
set PYTHONPATH=
```
**8.** Check appropriate version of vtk is installed and hdf5 is activated with 

```python
import vtk

print(vtk.VTK_VERSION)
try:
    vtk.vtkHDFReader
    vtk.vtkHDFWriter
    print("HDF5 support is available in your VTK build.")
except AttributeError:
    print("HDF5 support is NOT available in your VTK build.")
```
You can issue these lines of code on the terminal window by running python with command `python` or paste them into a `python.py` file and running the file with appropriate environment activated.

# Generating Sphinx documentation
From the root directory run the following commands [Powershell in windows]
- `./docs/make clean`
- `sphinx-apidoc -o ./docs/source/ ./fame`
- `./docs/make html`
- `./docs/make latexpdf` to generate pdf. 
<span style="color:red;">The appropriate latex compiler must be installed and available.</span>.

 Or bash in linux
- `make -C ./docs clean`
- `sphinx-apidoc -o ./docs/source/ ./fame`
- `make -C ./docs html`
- `make -C ./docs latexpdf`
<span style="color:red;">The appropriate latex compiler must be installed and available.</span>.

# Generating Docstring using VSCode extension

If you are using VS code, the [Python Docstring](https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring) extension can be used to auto-generate a docstring snippet once a function/class has been written. If you want the extension to generate docstrings in `Sphinx` format, you must set the `"autoDocstring.docstringFormat": "sphinx"` setting, under File > Preferences > Settings.

Note that it is best to write the docstrings once you have fully defined the function/class, as then the extension will generate the full dosctring. If you make any changes to the code once a docstring is generated, you will have to manually go and update the affected docstrings.

# FAME: FVM Based Laser Powder Bed Fusion Additive Manufacturing Process Simulation

## restructred text live preview on vscode

- need the extension `pip install esbonio` and 'esbonio'. After that make sure python path is manually setup if esbonio is having difficulty finding the python interpreter. You can do that by going to `File->Preference->Settings` and finding `Esbonio>Server:Python Path` For example: `C:\\msys64\\ucrt64\\bin\\python.exe` 

## restructured text cheatsheet

Details of rST is found [here](https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html)

## building the package

`python -m build`

## uploading the package to the test package host using "twine"
`twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose`
