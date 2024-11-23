import vtk

# Create a structured grid
dimensions = [10, 10, 1]  # Grid dimensions (nx, ny, nz)
points = vtk.vtkPoints()

for z in range(dimensions[2]):
    for y in range(dimensions[1]):
        for x in range(dimensions[0]):
            points.InsertNextPoint(x, y, z)

# Create a structured grid
structured_grid = vtk.vtkStructuredGrid()
structured_grid.SetDimensions(*dimensions)
structured_grid.SetPoints(points)

# Create scalar data
scalars = vtk.vtkFloatArray()
scalars.SetName("ScalarField")

# Assign scalar values (e.g., x^2 + y^2)
for z in range(dimensions[2]):
    for y in range(dimensions[1]):
        for x in range(dimensions[0]):
            scalar_value = x**2 + y**2
            scalars.InsertNextValue(scalar_value)

# Add the scalar data to the structured grid
structured_grid.GetPointData().SetScalars(scalars)

# Write the structured grid to a VTK file
# writer = vtk.vtkStructuredGridWriter()
# writer.SetFileName("structured_grid.vtk")
# writer.SetInputData(structured_grid)
# writer.Write()

# print("VTK file 'structured_grid.vtk' written successfully.")
