import vtk

# Define the points of the polygon
points = vtk.vtkPoints()
points.InsertNextPoint(1, 0, 0)
points.InsertNextPoint(1, 1, 0)
points.InsertNextPoint(1, 1, 1)
points.InsertNextPoint(1, 0, 1)

# Create a vtkPolygon and set its points
polygon = vtk.vtkPolygon()
polygon.GetPoints().DeepCopy(points)
polygon.GetPointIds().SetNumberOfIds(4)
for i in range(4):
    polygon.GetPointIds().SetId(i, i)

# Compute the area of the polygon
area = polygon.ComputeArea()
print(f"Computed Area: {area}")