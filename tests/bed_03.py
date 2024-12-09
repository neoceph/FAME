import vtk


import vtk
import numpy as np
from tqdm import tqdm
import h5py
import numpy as np


class StructuredMesh(vtk.vtkStructuredGrid):
    def __init__(self, bounds, divisions):
        """
        Initializes the StructuredMesh.
        
        Args:
            bounds (tuple): Bounds of the grid as ((x_min, x_max), (y_min, y_max), (z_min, z_max)).
            divisions (tuple): Number of divisions along x, y, z as (nx, ny, nz).
        """
        super().__init__()
        self.sharedCells = []
        self.cellCenters = []

        # Generate grid points
        self._generateGrid(bounds, divisions)

        # Compute cell centers and neighboring information
        self._computeCellCenter()
        self._computeNeighbors()

    def _generateGrid(self, bounds, divisions):
        """
        Generates the structured grid points and sets dimensions.
        """
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = bounds
        nx, ny, nz = divisions

        # Create points
        points = vtk.vtkPoints()
        dx = (x_max - x_min) / (nx - 1)
        dy = (y_max - y_min) / (ny - 1)
        dz = (z_max - z_min) / (nz - 1)

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    x = x_min + i * dx
                    y = y_min + j * dy
                    z = z_min + k * dz
                    points.InsertNextPoint(x, y, z)

        self.SetDimensions(nx, ny, nz)
        self.SetPoints(points)

    def _computeCellCenter(self):
        """
        Computes the centers of all cells using vtkCellCenters.
        """
        cell_centers_filter = vtk.vtkCellCenters()
        cell_centers_filter.SetInputData(self)
        cell_centers_filter.Update()

        cell_centers_output = cell_centers_filter.GetOutput()
        self.cell_centers = [
            cell_centers_output.GetPoint(i)
            for i in range(cell_centers_output.GetNumberOfPoints())
        ]

    def _computeNeighbors(self):
        """
        Computes shared cell information for all cells.
        """
        for cell_id in range(self.GetNumberOfCells()):
            currentCell = self.GetCell(cell_id)
            points = currentCell.GetPoints()

            vertices = [points.GetPoint(i) for i in range(points.GetNumberOfPoints())]
            shared_cells = []
            shared_vertices = []

            for other_cell_id in range(self.GetNumberOfCells()):
                if other_cell_id == cell_id:
                    continue

                other_cell = self.GetCell(other_cell_id)
                other_points = other_cell.GetPoints()
                shared_points = []

                for i in range(other_points.GetNumberOfPoints()):
                    if tuple(other_points.GetPoint(i)) in map(tuple, vertices):
                        shared_points.append(other_points.GetPoint(i))

                if len(shared_points) > 2:  # More than two vertices shared
                    shared_cells.append(other_cell_id)
                    shared_vertices.append(shared_points)

            self.sharedCells.append({
                "cell_id": cell_id,
                "shared_cells": shared_cells,
                "shared_vertices": shared_vertices
            })

    def getCellCenter(self, cell_id):
        """
        Retrieve the center of a specific cell.
        
        Args:
            cell_id (int): ID of the cell.
        
        Returns:
            tuple: Center of the cell (x, y, z).
        """
        if 0 <= cell_id < len(self.cell_centers):
            return self.cell_centers[cell_id]
        else:
            raise ValueError(f"Cell ID {cell_id} is out of range.")

    def getSharedCellsInfo(self, cell_id):
        """
        Retrieve shared cells information for a specific cell.
        
        Args:
            cell_id (int): ID of the cell.
        
        Returns:
            dict: Information about shared cells and shared vertices.
        """
        if 0 <= cell_id < len(self.sharedCells):
            return self.sharedCells[cell_id]
        else:
            raise ValueError(f"Cell ID {cell_id} is out of range.")



# Example Usage
if __name__ == "__main__":
    # Define bounds and divisions
    bounds = ((0, 10), (0, 5), (0, 2))  # ((x_min, x_max), (y_min, y_max), (z_min, z_max))
    divisions = (11, 6, 3)  # Number of divisions along x, y, z

    # Create a structured mesh
    mesh = StructuredMesh(bounds, divisions)

    # # Access cell centers and shared cells info
    # for cell_id in range(mesh.GetNumberOfCells()):
    #     print(f"Cell ID: {cell_id}")
    #     print(f"Cell Center: {mesh.get_cell_center(cell_id)}")
    #     print(f"Shared Cells Info: {mesh.get_shared_cells_info(cell_id)}")

    cell_id = 0
    for shared_cell, vertices in zip(mesh.getSharedCellsInfo(cell_id)["shared_cells"], mesh.getSharedCellsInfo(cell_id)["shared_vertices"]):
        print(f"Shared Cell ID: {shared_cell}")
        print(f"Shared Vertices: {vertices}")