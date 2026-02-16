import xml.etree.ElementTree as ET

# Check the benchmark file
file1 = 'benchmark_data/model_runs/runS0_specified_pressure_[0, 68.569]_final_output_Elements.vtu'
print(f"Checking: {file1}")
tree1 = ET.parse(file1)
root1 = tree1.getroot()

piece = root1.find('.//{http://www.vtk.org/vtkxml/version0.1}Piece')
if piece is not None:
    pd = piece.find('{http://www.vtk.org/vtkxml/version0.1}PointData')
    cd = piece.find('{http://www.vtk.org/vtkxml/version0.1}CellData')
    
    print(f"  PointData arrays: {len(pd) if pd is not None else 0}")
    if pd is not None:
        for arr in pd:
            name = arr.get('Name', 'unnamed')
            print(f"    - {name}")
    
    print(f"  CellData arrays: {len(cd) if cd is not None else 0}")
    if cd is not None:
        for arr in cd:
            name = arr.get('Name', 'unnamed')
            print(f"    - {name}")

print("\n" + "="*60)

# Now check the newer file
file2 = 'model_output/salt_wedge_benchmark/vtk_files/runS0_specified_pressure_[0, 68.569]_final_output.vtu'
print(f"Checking: {file2}")
tree2 = ET.parse(file2)
root2 = tree2.getroot()

piece2 = root2.find('.//{http://www.vtk.org/vtkxml/version0.1}Piece')
if piece2 is not None:
    pd2 = piece2.find('{http://www.vtk.org/vtkxml/version0.1}PointData')
    cd2 = piece2.find('{http://www.vtk.org/vtkxml/version0.1}CellData')
    
    print(f"  PointData arrays: {len(pd2) if pd2 is not None else 0}")
    if pd2 is not None:
        for arr in pd2:
            name = arr.get('Name', 'unnamed')
            print(f"    - {name}")
    
    print(f"  CellData arrays: {len(cd2) if cd2 is not None else 0}")
    if cd2 is not None:
        for arr in cd2:
            name = arr.get('Name', 'unnamed')
            print(f"    - {name}")
