import os
import vtk
import numpy as np
from .mesh import StructuredMesh

class MeshWriter:
    def __init__(self, mesh):
        """
        Initialize with a StructuredMesh instance.

        Args:
            mesh (StructuredMesh): The structured mesh object containing the grid and scalar data.
        """
        if not isinstance(mesh, StructuredMesh):
            raise TypeError("The provided mesh must be an instance of StructuredMesh.")
        self.mesh = mesh

    def _writeSingleVTS(self, output_file, variables):
        """
        Writes variables (scalar, vector, tensor) of the StructuredMesh to a .vts file.

        Args:
            output_file (str): Path to the output VTS file.
            variables (dict): Dictionary where keys are variable names and values are numpy arrays.
                              Each array should have the correct shape for the mesh dimensions.
        """
        if not output_file.endswith('.vts'):
            output_file += '.vts'

        if not hasattr(self.mesh, 'dimensions') or not hasattr(self.mesh, 'GetPoints'):
            raise ValueError("The provided mesh must have 'dimensions' and 'GetPoints' attributes.")

        nx, ny, nz = self.mesh.dimensions
        points = self.mesh.GetPoints()
        if points.GetNumberOfPoints() != nx * ny * nz:
            raise ValueError("Mismatch between vtkPoints size and grid dimensions.")

        for var_name, var_data in variables.items():
            if var_data.ndim == 1:
                if var_data.size != nx * ny * nz:
                    raise ValueError(f"Mismatch between scalar field '{var_name}' size and grid dimensions.")
                var_array = vtk.vtkDoubleArray()
                var_array.SetName(var_name)
                for value in var_data:
                    var_array.InsertNextValue(value)
                self.mesh.GetPointData().AddArray(var_array)
            elif var_data.ndim >= 2:
                if var_data.shape[0] != nx * ny * nz:
                    raise ValueError(f"Mismatch between field '{var_name}' size and grid dimensions.")
                if var_data.shape[1] == 2:
                    var_array = vtk.vtkDoubleArray()
                    var_array.SetName(var_name)
                    var_array.SetNumberOfComponents(2)
                    for value in var_data:
                        var_array.InsertNextTuple(value)
                    self.mesh.GetPointData().AddArray(var_array)
                elif var_data.shape[1] == 3:
                    var_array = vtk.vtkDoubleArray()
                    var_array.SetName(var_name)
                    var_array.SetNumberOfComponents(3)
                    for value in var_data:
                        var_array.InsertNextTuple(value)
                    self.mesh.GetPointData().AddArray(var_array)
                elif var_data.shape[1] == 4:
                    var_array = vtk.vtkDoubleArray()
                    var_array.SetName(var_name)
                    var_array.SetNumberOfComponents(4)
                    for value in var_data:
                        var_array.InsertNextTuple(value)
                    self.mesh.GetPointData().AddArray(var_array)
                elif var_data.shape[1] == 9:
                    var_array = vtk.vtkDoubleArray()
                    var_array.SetName(var_name)
                    var_array.SetNumberOfComponents(9)
                    for value in var_data:
                        var_array.InsertNextTuple(value)
                    self.mesh.GetPointData().AddArray(var_array)
                else:
                    raise ValueError(f"Unsupported number of components for field '{var_name}': {var_data.shape[1]}")
            else:
                raise ValueError(f"Unsupported data dimensions for variable '{var_name}': {var_data.ndim}")

        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(self.mesh)
        writer.Write()
        print(f"Structured mesh with variables written to {output_file}")

    def writeVTS(self, output_dir, variables, time=None, step=None):
        """
        Writes a single timestep data to a .vts file and updates the PVD file. For steady-state, it writes a single time step.

        Args:
            output_dir (str): Directory to save the output .vts file and PVD file.
            variables (dict): Dictionary of variables for the timestep.
            time (float, optional): Time value for the current timestep. Defaults to 0.0 for steady-state.
            step (int, optional): Step index for naming the .vts file. Defaults to 0 for steady-state.
        """
        os.makedirs(output_dir, exist_ok=True)
        pvd_file = os.path.join(output_dir, os.path.basename(output_dir) + '.pvd')

        if not os.path.exists(pvd_file):
            # Initialize a new PVD file if it doesn't exist
            with open(pvd_file, 'w', newline='') as f:
                f.write('<VTKFile type="Collection" version="0.1">\n')
                f.write('  <Collection>\n')
                f.write('  </Collection>\n')
                f.write('</VTKFile>\n')

        time = 0.0 if time is None else time
        step = 0 if step is None else step

        vts_file = os.path.join(output_dir, f"output_{step:04d}.vts")
        self._writeSingleVTS(vts_file, variables)

        # Update the PVD file
        with open(pvd_file, 'r+', newline='') as f:
            lines = f.readlines()
            insert_index = len(lines) - 2
            lines.insert(insert_index, f'    <DataSet timestep="{time}" group="" part="0" file="{os.path.basename(vts_file)}"/>\n')
            f.seek(0)
            f.writelines(lines)

        print(f"Updated PVD file: {pvd_file} with timestep {time} and file {vts_file}")
