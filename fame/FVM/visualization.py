import os
import vtk
import numpy as np
from .mesh import StructuredMesh, StructuredMesh1D

class MeshWriter:
    def __new__(cls, mesh):
        if isinstance(mesh, StructuredMesh1D):
            return super().__new__(MeshWriter1D)
        return super().__new__(MeshWriter3D)
    
    def __init__(self, mesh):
        """
        Initialize with a StructuredMesh instance.

        Args:
            mesh (StructuredMesh): The structured mesh object containing the grid and scalar data.
        """
        if not isinstance(mesh, StructuredMesh):
            raise TypeError("The provided mesh must be an instance of StructuredMesh.")
        self.mesh = mesh


class MeshWriter3D(MeshWriter):
    """
    Definition of 3D mesh writer.
    """
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

        points = self.mesh.GetPoints()
        num_points = points.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        if not points or num_points == 0:
            raise ValueError("The provided mesh must have valid vtkPoints.")

        for var_name, var_data in variables.items():
            var_array = vtk.vtkDoubleArray()
            var_array.SetName(var_name)

            if var_data.shape[0] == num_points:
                # Write to point data
                target = self.mesh.GetPointData()
            elif var_data.shape[0] == num_cells:
                # Write to cell data
                target = self.mesh.GetCellData()
            else:
                raise ValueError(
                    f"Mismatch between '{var_name}' size and grid dimensions. "
                    f"Expected {num_points} for points or {num_cells} for cells."
                )

            # Set number of components for scalar, vector, or tensor data
            if var_data.ndim == 1:
                var_array.SetNumberOfComponents(1)  # Explicit for scalar
                for value in var_data:
                    var_array.InsertNextValue(value)
            else:
                var_array.SetNumberOfComponents(var_data.shape[1])  # Vector or tensor
                for value in var_data:
                    var_array.InsertNextTuple(value)

            target.AddArray(var_array)

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


class MeshWriter1D(MeshWriter):
    def _writeSingleVTP(self, output_file, variables):
        if not output_file.endswith('.vtp'):
            output_file += '.vtp'

        points = self.mesh.GetPoints()
        num_points = points.GetNumberOfPoints()

        if not points or num_points == 0:
            raise ValueError("The provided mesh must have valid vtkPoints.")

        for var_name, var_data in variables.items():
            var_array = vtk.vtkDoubleArray()
            var_array.SetName(var_name)

            # Determine if data is for points or cells
            if var_data.shape[0] == num_points:
                target = self.mesh.GetPointData()
            elif var_data.shape[0] == self.mesh.GetNumberOfCells():
                target = self.mesh.GetCellData()
            else:
                raise ValueError(
                    f"Mismatch between '{var_name}' size and mesh. "
                    f"Expected {num_points} for points or {self.mesh.GetNumberOfCells()} for cells."
                )

            # Set components based on data shape (scalar, vector, tensor)
            if var_data.ndim == 1:
                var_array.SetNumberOfComponents(1)
                for value in var_data:
                    var_array.InsertNextValue(value)
            else:
                var_array.SetNumberOfComponents(var_data.shape[1])
                for value in var_data:
                    var_array.InsertNextTuple(value)

            target.AddArray(var_array)

        # Write the PolyData mesh to .vtp
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(self.mesh)
        writer.Write()

        print(f"PolyData mesh written to {output_file}")

    def writeVTS(self, output_dir, variables, time=None, step=None):
        os.makedirs(output_dir, exist_ok=True)
        pvd_file = os.path.join(output_dir, os.path.basename(output_dir) + '.pvd')

        # Create new PVD file if it doesn't exist
        if not os.path.exists(pvd_file):
            with open(pvd_file, 'w') as f:
                f.write('<VTKFile type="Collection" version="0.1">\n')
                f.write('  <Collection>\n')
                f.write('  </Collection>\n')
                f.write('</VTKFile>\n')

        time = 0.0 if time is None else time
        step = 0 if step is None else step

        vtp_file = os.path.join(output_dir, f"output_{step:04d}.vtp")
        self._writeSingleVTP(vtp_file, variables)

        # Update the PVD file with the new VTP entry
        with open(pvd_file, 'r+') as f:
            lines = f.readlines()
            insert_index = len(lines) - 2
            lines.insert(insert_index, f'    <DataSet timestep="{time}" file="{os.path.basename(vtp_file)}"/>\n')
            f.seek(0)
            f.writelines(lines)

        print(f"Updated PVD file: {pvd_file}")
