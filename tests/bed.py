## vtk methods and classes to watch
# vtk.StructuredGridConnectivity
# vtk.PStructedGridConnectivity

import vtk

# Create a sample structured grid
structuredGrid = vtk.vtkStructuredGrid()
points = vtk.vtkPoints()

# Define grid dimensions (e.g., 2x2x2 structured grid)
nx, ny, nz = 3, 3, 3
structuredGrid.SetDimensions(nx, ny, nz)

# Add points to the structured grid
for z in range(nz):
    for y in range(ny):
        for x in range(nx):
            points.InsertNextPoint(x, y, z)

structuredGrid.SetPoints(points)

# Initialize arrays to store face areas and normals
face_areas = vtk.vtkDoubleArray()
face_areas.SetName("Face Areas")

face_normals = vtk.vtkDoubleArray()
face_normals.SetName("Face Normals")
face_normals.SetNumberOfComponents(3)  # Normals are 3D vectors

# Temporary variables for calculation
cell = vtk.vtkGenericCell()
polygon = vtk.vtkPolygon()

for i in range(structuredGrid.GetNumberOfCells()):
    structuredGrid.GetCell(i, cell)
    
    for face_id in range(cell.GetNumberOfFaces()):
        face = cell.GetFace(face_id)
        if face is not None:
            # Extract points and IDs from the face
            points_of_face = face.GetPoints()
            num_points = points_of_face.GetNumberOfPoints()

            if num_points < 3:
                print(f"Skipping degenerate face with fewer than 3 points: {num_points}")
                face_areas.InsertNextValue(0.0)
                face_normals.InsertNextTuple([0.0, 0.0, 0.0])
                continue

            point_ids = [face.GetPointId(j) for j in range(num_points)]
            normal = [0.0, 0.0, 0.0]

            # Debug: Print points and IDs
            print(f"Processing face {face_id} of cell {i}")
            print(f"Number of points: {num_points}")
            for j in range(num_points):
                print(f"Point {j}: {points_of_face.GetPoint(j)}")

            # Compute the area and normal of the face
            try:
                area = vtk.vtkPolygon.ComputeArea(
                    points_of_face,
                    num_points,
                    point_ids,
                    normal
                )
                if area != area:  # Check for nan
                    print("Computed area is nan. Skipping this face.")
                    area = 0.0
                    normal = [0.0, 0.0, 0.0]

                face_areas.InsertNextValue(area)
                face_normals.InsertNextTuple(normal)
            except Exception as e:
                print(f"Error computing area or normal for face {face_id}: {e}")
                face_areas.InsertNextValue(0.0)
                face_normals.InsertNextTuple([0.0, 0.0, 0.0])

# Attach the computed face areas and normals to the structured grid for visualization
structuredGrid.GetCellData().AddArray(face_areas)
structuredGrid.GetCellData().AddArray(face_normals)

# Write the structured grid to a file for visualization (optional)
writer = vtk.vtkXMLStructuredGridWriter()
writer.SetFileName("structured_grid_with_face_areas_and_normals.vts")
writer.SetInputData(structuredGrid)
writer.Write()

print("Face areas and normals have been computed and added to the structured grid.")
