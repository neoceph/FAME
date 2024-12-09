import vtk

# Create a vtkStructuredGrid
structured_grid = vtk.vtkStructuredGrid()

# Define the dimensions of the grid
dims = [4, 4, 3]  # Grid dimensions (3x3x2)
structured_grid.SetDimensions(dims)

# Create points for the structured grid
points = vtk.vtkPoints()
for k in range(dims[2]):
    for j in range(dims[1]):
        for i in range(dims[0]):
            points.InsertNextPoint(i, j, k)

structured_grid.SetPoints(points)

# Compute the cell centers
cell_centers = vtk.vtkCellCenters()
cell_centers.SetInputData(structured_grid)
cell_centers.Update()

# Extract the centers
centers = cell_centers.GetOutput()

# Print the cell centers
num_centers = centers.GetNumberOfPoints()
print("Cell Centers:")
for i in range(num_centers):
    x, y, z = centers.GetPoint(i)
    print(f"Cell {i}: Center at ({x:.2f}, {y:.2f}, {z:.2f})")
