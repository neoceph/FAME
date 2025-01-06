import vtk
import numpy as np
import scipy.sparse as sp

from scipy.spatial import ConvexHull
from tqdm import tqdm


class StructuredMesh:
    def __new__(cls, bounds, divisions):
        """
        Dynamically instantiate the correct subclass based on the dimensionality.
        """
        if len(divisions) == 1:
            return vtk.vtkPolyData.__new__(StructuredMesh1D)
        else:
            return vtk.vtkStructuredGrid.__new__(StructuredMesh3D)    
    
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
        self.divisions = divisions
        # self.is_1D = len(divisions) == 1

        # Create vtkStructuredGrid or vtkPolyData depending on 1D or 3D mesh
        # self.mesh = vtk.vtkPolyData() if self.is_1D else vtk.vtkStructuredGrid()

        # Generate grid points
        self._generateGrid(bounds, divisions)
        self._computeCellFaces()
        self._computeCellCenter()
        self._computeNeighbors()
        
        self.numCells = self.GetNumberOfCells()

        self.A = sp.lil_matrix((self.numCells, self.numCells))  # Use LIL format for construction
        self.b = np.zeros(self.numCells)


class StructuredMesh3D(StructuredMesh, vtk.vtkStructuredGrid):
    
    def __init__(self, bounds, divisions):
        vtk.vtkStructuredGrid.__init__(self)
        super().__init__(bounds, divisions)

    def GetNumberOfCells(self):
        return super().GetNumberOfCells()
    
    def _generateGrid(self, bounds, divisions):
        """
        Generates the structured grid points and sets dimensions.
        """
        # Create points
        points = vtk.vtkPoints()
        
        (x_min, x_max), (y_min, y_max), (z_min, z_max) = bounds
        div_x, div_y, div_z = divisions

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
        self.cellCenters = [
            cell_centers_output.GetPoint(i)
            for i in range(cell_centers_output.GetNumberOfPoints())
        ]

    def _computeNeighbors(self):
        """
        Computes shared cell information for all cells.
        """
        self.sharedCells = []
        for cell_id in range(self.GetNumberOfCells()):
            currentCell = self.GetCell(cell_id)
            points = currentCell.GetPoints()

            vertices = {tuple(points.GetPoint(i)) for i in range(points.GetNumberOfPoints())}
            shared_cells = []
            shared_faces = set()

            for other_cell_id in range(self.GetNumberOfCells()):
                if other_cell_id == cell_id:
                    continue

                other_cell = self.GetCell(other_cell_id)
                other_points = other_cell.GetPoints()
                other_vertices = {tuple(other_points.GetPoint(i)) for i in range(other_points.GetNumberOfPoints())}
                
                shared_points = vertices.intersection(other_vertices)

                if len(shared_points) > 2:
                    shared_cells.append(other_cell_id)
                    shared_faces.update(face_id for face_id, face_points in self.faces.items()
                                        if set(self.GetPoint(pt) for pt in face_points) == shared_points)

            # Identify boundary faces by subtracting shared faces from all faces of the cell
            all_faces = {
                face_id for face_id, face_points in self.faces.items()
                if set(self.GetPoint(pt) for pt in face_points).issubset(vertices)
            }
            boundary_faces = all_faces - shared_faces

            self.sharedCells.append({
                "cell_id": cell_id,
                "shared_cells": shared_cells,
                "shared_faces": list(shared_faces),
                "boundary_faces": list(boundary_faces)
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
        if 0 <= cell_id < len(self.cellCenters):
            return self.cellCenters[cell_id]
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

    def getFaceByCenter(self, center, tolerance=1e-6):
        """
        Retrieve all face IDs sorted by proximity to a center point, optionally filtered by a tolerance.

        Args:
            center (tuple or list): Target center point (x, y, z).
            tolerance (float, optional): Distance tolerance. Defaults to 1e-9.

        Returns:
            list: List of face IDs sorted by proximity to the center point.
        """
        if not isinstance(center, (tuple, list)) or len(center) != 3:
            raise ValueError("Center must be a tuple or list of length 3.")
        
        target = np.array(center)

        distances = {fid: np.linalg.norm(np.array(self.faceCenters[fid]) - target)
                    for fid in self.faceCenters}

        return sorted([fid for fid, dist in distances.items() if dist <= tolerance], key=lambda fid: distances[fid])


    def getFacesByCoordinates(self, x=None, y=None, z=None, tolerance=None):
        """
        Retrieve face IDs by proximity to x, y, or z coordinates, optionally within a tolerance.

        Args:
            x (float, optional): x-coordinate to filter faces.
            y (float, optional): y-coordinate to filter faces.
            z (float, optional): z-coordinate to filter faces.
            tolerance (float, optional): Distance tolerance. If provided, returns all faces within the tolerance along specified axes.

        Returns:
            list: List of face IDs within the tolerance of the provided coordinates.
        
        Raises:
            ValueError: If none of x, y, z are provided.
        """
        if x is None and y is None and z is None:
            raise ValueError("At least one of x, y, or z must be provided.")
        
        # Default tolerance if not specified
        if tolerance is None:
            tolerance = 1e-6

        matching_faces = []
        
        for fid, center in self.faceCenters.items():
            match = True
            if x is not None and not (abs(center[0] - x) <= tolerance):
                match = False
            if y is not None and not (abs(center[1] - y) <= tolerance):
                match = False
            if z is not None and not (abs(center[2] - z) <= tolerance):
                match = False
            
            if match:
                matching_faces.append(fid)
        
        return matching_faces

    def getCellIdByFaceId(self, face_id):
        """
        Retrieve the cell ID that owns a specific face ID.

        Args:
            face_id (int): ID of the face to search for.

        Returns:
            int: The cell ID that owns the specified face.
        
        Raises:
            ValueError: If the face ID does not belong to any cell.
        """
        for cell in self.sharedCells:
            if face_id in cell['shared_faces'] or face_id in cell['boundary_faces']:
                return cell['cell_id']
        
        raise ValueError(f"Face ID {face_id} does not belong to any cell.")


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

    def getCellVolume(self, cell_id):
        """
        Calculate the volume of a cell using the Gauss Divergence Theorem.

        Args:
            cell_id (int): ID of the cell.

        Returns:
            float: Volume of the cell.
        """
        if cell_id < 0 or cell_id >= self.GetNumberOfCells():
            raise ValueError(f"Cell ID {cell_id} is out of range.")
        
        # Retrieve shared and boundary faces for the cell
        shared_faces_info = self.sharedCells[cell_id]
        face_ids = shared_faces_info['shared_faces'] + shared_faces_info['boundary_faces']
        
        volume = 0.0
        for face_id in face_ids:
            face_pts = self.faces[face_id]

            # Calculate face area and normal using existing calculateArea method
            vtk_points = vtk.vtkPoints()
            for pt_id in face_pts:
                vtk_points.InsertNextPoint(self.GetPoint(pt_id))
            
            area = self.calculateArea(vtk_points)
            
            volume += area
        
        volume /= len(face_ids)

        return abs(volume)


class StructuredMesh1D(StructuredMesh, vtk.vtkPolyData):
    def __init__(self, bounds, divisions):
        vtk.vtkPolyData.__init__(self)
        super().__init__(bounds, divisions)

        # Define faceArea as a float variable specific to 1D mesh
        self.faceArea = 1.0

    def GetNumberOfCells(self):
        return self.GetNumberOfLines()

    def _generateGrid(self, bounds, divisions):
        """
        Generates the structured grid points and sets dimensions.
        """
        # Create points
        points = vtk.vtkPoints()
        
        # 1D mesh: Generate points along a line
        (x_min, x_max) = bounds
        div_x = divisions[0]
        dx = (x_max - x_min) / div_x

        for i in range(div_x + 1):
            x = x_min + i * dx
            points.InsertNextPoint(x, 0, 0)  # 1D is along x-axis, y/z set to 0

        # Create vtkPolyData for 1D mesh
        lines = vtk.vtkCellArray()
        for i in range(div_x):
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(line)

        self.SetPoints(points)
        self.SetLines(lines)
        
    def _computeCellFaces(self):

        # Directly use point IDs to represent unique faces
        self.faces = {}  # Reset face storage
        self.faceCenters = {}

        num_points = self.GetNumberOfPoints()

        for i in range(num_points):
            # Use the point ID as the face ID
            point_id = i
            self.faces[point_id] = (point_id,)
            self.faceCenters[point_id] = self.GetPoint(point_id)

    def _computeCellCenter(self):
        """
        Computes the centers of all cells using vtkCellCenters.
        """
        # For 1D, compute center of each line segment
        self.cellCenters = []
        for i in range(self.GetNumberOfCells()):
            line = self.GetCell(i)
            p1 = np.array(self.GetPoint(line.GetPointId(0)))
            p2 = np.array(self.GetPoint(line.GetPointId(1)))
            center = (p1 + p2) / 2
            self.cellCenters.append(center)

    def _computeNeighbors(self):
        """
        Computes shared cell information for all cells.
        """
        self.sharedCells = []
        num_cells = self.GetNumberOfCells()

        # Iterate over cells (line segments)
        for cell_id in range(num_cells):
            line = self.GetCell(cell_id)
            p1 = line.GetPointId(0)  # Start point of line segment
            p2 = line.GetPointId(1)  # End point of line segment

            # Shared cells and faces
            shared_cells = []
            shared_faces = set()

            if cell_id > 0:
                shared_cells.append(cell_id - 1)  # Left neighbor
                shared_faces.add(p1)  # Shared face at left endpoint

            if cell_id < num_cells - 1:
                shared_cells.append(cell_id + 1)  # Right neighbor
                shared_faces.add(p2)  # Shared face at right endpoint

            # Boundary faces
            boundary_faces = set()
            if cell_id == 0:  # Leftmost cell
                boundary_faces.add(p1)
            if cell_id == num_cells - 1:  # Rightmost cell
                boundary_faces.add(p2)

            # Store results
            self.sharedCells.append({
                "cell_id": cell_id,
                "shared_cells": shared_cells,
                "shared_faces": list(shared_faces),
                "boundary_faces": list(boundary_faces)
            })

    def _computeCellFaces(self):
        # Directly use point IDs to represent unique faces
        self.faces = {}  # Reset face storage
        self.faceCenters = {}

        num_points = self.GetNumberOfPoints()

        for i in range(num_points):
            # Use the point ID as the face ID
            point_id = i
            self.faces[point_id] = (point_id,)
            self.faceCenters[point_id] = self.GetPoint(point_id)

    def GetDimensions(self):
        """
        Returns the dimensions for the 1D mesh.
        
        For a 1D mesh, the x-dimension is derived from the number of cells,
        and y and z dimensions are always 1.
        
        Returns:
            tuple: (nx, ny, nz)
        """
        n = self.GetNumberOfCells() + 1  # Number of points is cells + 1
        return n

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
        
    def getFacesByCoordinates(self, x=None, tolerance=1e-6):
        """
        Retrieve face IDs by proximity to x-coordinate in 1D mesh.

        Args:
            x (float, optional): x-coordinate to filter faces.
            tolerance (float, optional): Distance tolerance for matching faces.

        Returns:
            list: List of face IDs within the tolerance of the provided x-coordinate.

        Raises:
            ValueError: If x is not provided.
        """
        if x is None:
            raise ValueError("x-coordinate must be provided for 1D mesh.")

        matching_faces = []
        for fid, center in self.faceCenters.items():
            if abs(center[0] - x) <= tolerance:
                matching_faces.append(fid)
        
        return matching_faces
        
    def calculateArea(self, vtk_points, includeNormal=False):
        """
        Return the constant face area if nothing is specified
        """

        return self.faceArea
    
    def getCellVolume(self, cell_id):
        """
        Calculate the volume of a 1D cell using its length and face area.

        Args:
            cell_id (int): ID of the cell.

        Returns:
            float: Volume of the cell.
        
        Raises:
            ValueError: If the cell ID is out of range.
        """
        if cell_id < 0 or cell_id >= self.GetNumberOfCells():
            raise ValueError(f"Cell ID {cell_id} is out of range.")
        
        # Get the line segment representing the cell
        line = self.GetCell(cell_id)
        p1 = np.array(self.GetPoint(line.GetPointId(0)))
        p2 = np.array(self.GetPoint(line.GetPointId(1)))

        # Calculate the length of the line segment
        cell_length = np.linalg.norm(p2 - p1)

        # Calculate volume (length * cross-sectional area)
        volume = cell_length * self.faceArea

        return volume