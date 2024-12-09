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

# Define a helper function to find neighboring cells that share a face
def find_shared_face_cells(grid, cell_id):
    dims = grid.GetDimensions()
    i, j, k = [cell_id % dims[0], (cell_id // dims[0]) % dims[1], cell_id // (dims[0] * dims[1])]
    neighbors = []

    directions = [
        (-1, 0, 0), (1, 0, 0),  # Left and Right
        (0, -1, 0), (0, 1, 0),  # Front and Back
        (0, 0, -1), (0, 0, 1)   # Bottom and Top
    ]

    for di, dj, dk in directions:
        ni, nj, nk = i + di, j + dj, k + dk
        if 0 <= ni < dims[0] - 1 and 0 <= nj < dims[1] - 1 and 0 <= nk < dims[2] - 1:
            neighbor_id = ni + (nj * (dims[0] - 1)) + (nk * (dims[0] - 1) * (dims[1] - 1))
            neighbors.append(neighbor_id)
    return neighbors

# Iterate through cells and find neighbors
num_cells = structured_grid.GetNumberOfCells()
for cell_id in range(num_cells):
    neighbors = find_shared_face_cells(structured_grid, cell_id)
    print(f"Cell {cell_id} shares faces with cells: {neighbors}")
