import vtk
import numpy as np

from scipy.spatial import ConvexHull
from tqdm import tqdm


class StructuredMesh(vtk.vtkStructuredGrid):
    def __init__(self, bounds, divisions):
        """
        Initializes the StructuredMesh.
        
        Args:
            bounds (tuple): Bounds of the grid as ((x_min, x_max), (y_min, y_max), (z_min, z_max)).
            divisions (tuple): Number of divisions along x, y, z as (div_x, div_y, div_z).
        """
        super().__init__()
        self.sharedCells = []
        self.cellCenters = []
        self.faces = {}
        self.faceCenters = {}
    

        # Generate grid points
        self._generateGrid(bounds, divisions)

        # Compute cell centers, neighborcell, cell faces
        self._computeCellFaces()
        self._computeCellCenter()
        self._computeNeighbors()

    def _generateGrid(self, bounds, divisions):
        """
        Generates the structured grid points and sets dimensions.
        """
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = bounds
        div_x, div_y, div_z = divisions

        # Create points
        points = vtk.vtkPoints()
        dx = (x_max - x_min) / div_x
        dy = (y_max - y_min) / div_y
        dz = (z_max - z_min) / div_z

        nx = div_x + 1
        ny = div_y + 1
        nz = div_z + 1

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

            vertices = {tuple(points.GetPoint(i)) for i in range(points.GetNumberOfPoints())}
            shared_cells = []
            shared_faces = []

            for other_cell_id in range(self.GetNumberOfCells()):
                if other_cell_id == cell_id:
                    continue

                other_cell = self.GetCell(other_cell_id)
                other_points = other_cell.GetPoints()
                other_vertices = {tuple(other_points.GetPoint(i)) for i in range(other_points.GetNumberOfPoints())}
                
                shared_points = vertices.intersection(other_vertices)

                if len(shared_points) > 2:
                    shared_cells.append(other_cell_id)
                    shared_faces.extend(face_id for face_id, face_points in self.faces.items()
                                        if set(self.GetPoint(pt) for pt in face_points) == shared_points)

            self.sharedCells.append({
                "cell_id": cell_id,
                "shared_cells": shared_cells,
                "shared_faces": shared_faces
            })

    def _computeCellFaces(self):
        face_set = set()
        face_id = 0

        cell_iter = self.NewCellIterator()

        while not cell_iter.IsDoneWithTraversal():
            cell = vtk.vtkGenericCell()
            cell_iter.GetCell(cell)
            faces = self._extractCellFaces(cell, hexahedron=True)

            for face in faces:
                face_tuple = tuple(sorted(face))
                if face_tuple not in face_set:
                    face_set.add(face_tuple)
                    self.faces[face_id] = face_tuple
                    center = np.mean([self.GetPoint(idx) for idx in face_tuple], axis=0)
                    self.faceCenters[face_id] = tuple(center)
                    face_id += 1

            cell_iter.GoToNextCell()

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

    def getFaceById(self, face_id):
        """
        Retrieve face points by face ID.
        """
        return self.faces.get(face_id, None)

    def getFaceByCenter(self, center):
        """
        Retrieve face ID by proximity to a center point.
        """
        center = np.array(center)
        return min(self.faceCenters, key=lambda fid: np.linalg.norm(np.array(self.faceCenters[fid]) - center))

    def listFacesByPoint(self, point_id):
        """
        Retrieve all faces associated with a given point.
        """
        return self.pointFaces.get(point_id, [])        

    def _extractCellFaces(self, cell, hexahedron=False):
        """
        Extracts the six faces from a hexahedral cell.
        Returns a list of 4-point faces using global point IDs from vtkStructuredGrid.
        Hexahedral Cell Point IDs (VTK Order):
          7-------6
         /|      /|
        4-------5 |
        | |     | |
        | 3-----|-2
        |/      |/
        0-------1        
        """
        point_ids = cell.GetPointIds()
        
        hexahedron_faces = [
            [0, 1, 5, 4],   # -Y Face
            [1, 2, 6, 5],   # +X Face
            [2, 3, 7, 6],   # +Y Face
            [3, 0, 4, 7],   # -X Face
            [0, 1, 2, 3],   # -Z Face
            [4, 5, 6, 7]    # +Z Face
        ]

        faces = hexahedron_faces
        face_points = []
        for face in faces:
            face_points.append([point_ids.GetId(i) for i in face])

        return face_points

    def calculateArea(self, vtk_points, includeNormal=False):
        """
        Calculate the area of a planar polygon using the Shoelace Formula.

        Args:
            vtk_points (vtk.vtkPoints): vtkPoints object containing 3D coordinates of the polygon vertices.

        Returns:
            float: The computed area of the planar polygon.
        """
        # Convert vtkPoints to a numpy array
        num_points = vtk_points.GetNumberOfPoints()
        if num_points < 3:
            raise ValueError("At least 3 points are required to calculate the area.")

        points = np.array([vtk_points.GetPoint(i) for i in range(num_points)])

        # Step 1: Verify if points are planar
        # Compute the normal vector of the plane using the first three points
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        normal = np.cross(v1, v2)
        if np.linalg.norm(normal) == 0:
            raise ValueError("Points are collinear and do not form a polygon.")
        normal = normal / np.linalg.norm(normal)  # Normalize the normal vector

        # Check planarity by ensuring all points lie on the same plane
        for point in points[3:]:
            if not np.isclose(np.dot(point - points[0], normal), 0, atol=1e-6):
                raise ValueError("Points are not planar and do not form a valid polygon.")

        # Step 2: Project points onto the dominant 2D plane
        # Drop the axis with the largest normal component
        dominant_axis = np.argmax(np.abs(normal))  # Choose axis to drop
        projected_points = np.delete(points, dominant_axis, axis=1)

        # Ensure projected points are ordered counterclockwise
        hull = ConvexHull(projected_points)  # Ensures proper ordering
        ordered_points = projected_points[hull.vertices]

        # Step 3: Calculate the area using the Shoelace Formula
        x = ordered_points[:, 0]
        y = ordered_points[:, 1]
        area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

        if includeNormal:
            return area, normal
        else:
            return area
