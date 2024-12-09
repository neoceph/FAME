import vtk
import time

def create_structured_grid(dims):
    """Create a vtkStructuredGrid with specified dimensions."""
    structured_grid = vtk.vtkStructuredGrid()
    structured_grid.SetDimensions(dims)
    
    # Generate points for the grid
    points = vtk.vtkPoints()
    for z in range(dims[2]):
        for y in range(dims[1]):
            for x in range(dims[0]):
                points.InsertNextPoint(x, y, z)
    
    structured_grid.SetPoints(points)
    return structured_grid

def find_shared_vertices(grid):
    """Identify cells that share vertices and explicitly list shared vertices."""
    shared_info = {}  # Dictionary to store shared information
    num_cells = grid.GetNumberOfCells()
    point_to_cells = {}  # Map points to the cells they belong to

    # Build a point-to-cell map
    for cell_id in range(num_cells):
        cell = grid.GetCell(cell_id)
        for i in range(cell.GetNumberOfPoints()):
            point_id = cell.GetPointId(i)
            if point_id not in point_to_cells:
                point_to_cells[point_id] = []
            point_to_cells[point_id].append(cell_id)

    # Find shared vertices for each pair of cells
    for cell_id in range(num_cells):
        current_cell = grid.GetCell(cell_id)
        for i in range(current_cell.GetNumberOfPoints()):
            point_id = current_cell.GetPointId(i)
            for neighbor_cell_id in point_to_cells[point_id]:
                if neighbor_cell_id != cell_id:  # Skip self
                    shared_vertices = shared_info.setdefault((cell_id, neighbor_cell_id), set())
                    shared_vertices.add(point_id)

    return shared_info

def print_shared_info(shared_info):
    """Print shared vertices information for each pair of cells."""
    detailed_shared_info = {}

    # Organize shared info for better printing
    for cell_pair, vertices in shared_info.items():
        cell1, cell2 = cell_pair
        if cell1 not in detailed_shared_info:
            detailed_shared_info[cell1] = {}
        detailed_shared_info[cell1][cell2] = vertices

    # Print detailed shared information
    for cell1, neighbors in detailed_shared_info.items():
        for cell2, vertices in neighbors.items():
            vertices_str = ", ".join(map(str, sorted(vertices)))
            print(f"Cell {cell1} shares vertices {vertices_str} with cell {cell2}")

# Main Execution
if __name__ == "__main__":
    # Define grid dimensions (number of points in x, y, z directions)
    dimensions = [4, 4, 3]  # A 4x4x4 grid of cells

    # Measure grid creation time
    start_time = time.time()
    structured_grid = create_structured_grid(dimensions)
    grid_creation_time = time.time() - start_time

    # Measure shared vertices computation time
    start_time = time.time()
    shared_info = find_shared_vertices(structured_grid)
    computation_time = time.time() - start_time

    # Measure result printing time
    start_time = time.time()
    # print_shared_info(shared_info)
    printing_time = time.time() - start_time

    # Print timing results
    print("\nTiming Results:")
    print(f"Grid Creation Time: {grid_creation_time:.6f} seconds")
    print(f"Shared Vertices Computation Time: {computation_time:.6f} seconds")
    print(f"Result Printing Time: {printing_time:.6f} seconds")
