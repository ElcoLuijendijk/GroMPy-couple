# ParaView Visualization Guide for GroMPy-couple VTK Output

A comprehensive guide to visualizing GroMPy-couple simulation results in ParaView.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Opening VTK Files](#opening-vtk-files)
4. [Variable Visualization](#variable-visualization)
5. [Creating Visualizations](#creating-visualizations)
6. [Advanced Features](#advanced-features)
7. [Troubleshooting](#troubleshooting)

---

## Installation

### ParaView

ParaView is a free, open-source visualization application maintained by Kitware.

**Download and Install:**
- Visit https://www.paraview.org/download/
- Choose the appropriate installer for your OS (Windows, macOS, Linux)
- Follow the installation wizard
- No additional dependencies needed

**Minimum Requirements:**
- ParaView version 5.4+ (any recent version will work)
- 2 GB RAM recommended (4 GB for large models)
- OpenGL 2.1+ support

### Python VTK Support (for GroMPy-couple)

Ensure GroMPy-couple can generate VTK files:

```bash
# Install native VTK bindings (recommended)
pip install vtkmodules

# Or use meshio as fallback
pip install meshio
```

---

## Quick Start

### 1. Generate VTK Files

In your model parameters file (`model_input/model_parameters.py`):

```python
class ModelOptions:
    save_vtk_files = True
    save_vtk_files_all_steps = False  # Set True for per-timestep output
```

Run your model:

```bash
python grompy.py model_input/model_parameters.py
```

VTK files will be saved to `model_output/<scenario>/vtk_files/`

### 2. Open in ParaView

```bash
# Open all VTK files in a directory
paraview model_output/*/vtk_files/

# Or open specific file
paraview model_output/terrestrial_steady/vtk_files/runS0_final_output_Elements.vtu
```

ParaView will launch with your visualization ready.

---

## Opening VTK Files

### Method 1: Command Line

```bash
# Single file
paraview model_output/scenario/vtk_files/run_final_output_Elements.vtu

# All files in directory (loads as file series)
paraview model_output/scenario/vtk_files/*.vtu

# With specific variable visualization
paraview --data=model_output/scenario/vtk_files/*.vtu
```

### Method 2: ParaView GUI

1. **Open ParaView**
2. **File → Open**
3. Navigate to `model_output/<scenario>/vtk_files/`
4. Select desired `.vtu` file
5. Click "Open"
6. In the "Open Data With..." dialog, select "XML Unstructured Grid Reader"
7. Click "OK"

### Method 3: Drag and Drop

- Drag VTK files directly into ParaView window
- ParaView will automatically detect the format

---

## Variable Visualization

### Available Variables

GroMPy-couple outputs 18 variables per simulation:

#### Primary Fields (Most Useful for Visualization)

| Variable | Units | Purpose | Notes |
|----------|-------|---------|-------|
| **pressure** | Pa | Fluid pressure | Color code for pressure distribution |
| **concentration** | ppt | Solute concentration (salinity) | Visualize salt distribution |
| **hydraulic_head** | m | Hydraulic head | Alternative to pressure |
| **velocity_x, velocity_y** | m/s | Velocity components | Plot with arrows for flow direction |
| **flux_magnitude** | m/s | Velocity magnitude | Use for streamlines/glyphs |

#### Permeability & Flux

| Variable | Units | Purpose | Notes |
|----------|-------|---------|-------|
| **permeability_xx, permeability_yy** | m² | Permeability tensor (diagonal) | Visualize heterogeneity |
| **nodal_flux** | m³/s | Surface flux | Non-zero only on boundaries |

#### Boundary Masks (Boolean: 0/1)

| Variable | Purpose | Usage |
|----------|---------|-------|
| **surface_mask** | Top surface identification | Use with threshold filter |
| **sea_surface_mask** | Seawater surface (coastal models) | Coastal aquifer studies |
| **specified_pressure_bnd** | Pressure boundary locations | Boundary condition visualization |
| **active_seepage_bnd** | Seepage boundary locations | Seepage face identification |
| **recharge_bnd** | Recharge boundary locations | Recharge area visualization |
| **active_concentration_bnd** | Concentration boundary locations | Salt source visualization |

### Visualizing Pressure

1. **Load VTK file** (as described above)
2. In the left panel, check the variable list under "Properties"
3. Select **"pressure"** from the "Coloring" dropdown (top toolbar)
4. Adjust color scale:
   - Right-click on pressure variable
   - Select "Edit Color Map"
   - Choose color scheme (e.g., "Cool to Warm", "Viridis")
5. Use color bar legend to interpret values

### Visualizing Concentration (Salinity)

1. Select **"concentration"** from the "Coloring" dropdown
2. Use contrasting colors to distinguish fresh/salt water:
   - Blue for fresh water (0 ppt)
   - Red for salt water (35 ppt)
3. Useful for studying saltwater intrusion

### Visualizing Velocity Field

#### Method 1: Color by Magnitude

1. Select **"flux_magnitude"** from coloring dropdown
2. Adjust color scale as needed

#### Method 2: Glyphs (Arrows)

This shows flow direction:

1. **Filters → Visualization → Glyph**
2. Configure glyph properties:
   - **Glyph Type:** Arrow
   - **Vectors:** flux_magnitude (or velocity_x/y)
   - **Scale Mode:** Scale by Magnitude
   - **Scale Factor:** Adjust as needed (start with 10)
3. Click **Apply**

#### Method 3: Streamlines

Shows flow paths:

1. **Filters → Visualization → Stream Tracer**
2. Configure:
   - **Vectors:** flux_magnitude
   - **Seed Type:** High Resolution Line Source
   - **Line locations:** Define line across model
3. Click **Apply**
4. Adjust opacity, color, and width in Properties panel

---

## Creating Visualizations

### Single Variable Visualization

**Example: Pressure Distribution**

1. Open VTK file in ParaView
2. Select **"pressure"** from coloring dropdown
3. Adjust color scale (Right-click color bar → Edit)
4. Optional: Add **Surface Representation**
   - Top toolbar: Change "Surface" dropdown to "Wireframe" for mesh visibility
   - Or use "Surface with Edges" for both

### Multiple Variables (Simultaneous Display)

**Example: Pressure + Flow Velocity**

1. Load VTK file
2. **Duplicate representation** of the data:
   - Right-click data in Pipeline Browser (left panel)
   - Select "Duplicate"
3. On first representation: Color by **"pressure"**
4. On second representation: Apply **Glyph** filter with arrows
5. Adjust opacity of second representation so both are visible

### Threshold Visualization

**Example: Show only seawater regions**

1. Load VTK file
2. **Filters → Threshold**
3. Configure:
   - **Column:** sea_surface_mask
   - **Value Range:** 0.5 to 1.5 (selects cells with mask=1)
4. Click **Apply**
5. Color by **"concentration"** to see salinity distribution in seawater only

### Slice Visualization

**Cross-section through the model**

1. Load VTK file
2. **Filters → Slice**
3. Configure:
   - **Plane Location:** Set Y or X coordinate for cross-section
   - **Plane Normal:** Choose cutting direction
4. Click **Apply**
5. This creates a 2D surface through the 3D domain
6. Color by desired variable

---

## Advanced Features

### Exporting High-Quality Images

**Publication-ready figures:**

1. Prepare visualization as desired
2. **File → Export Scene**
3. Choose format:
   - **PNG** (recommended, good quality)
   - **PDF** (vector format, very clean)
   - **EPS** (scientific publishing)
4. Configure resolution:
   - Set higher pixel dimensions for printing (e.g., 3000×2000)
5. Adjust settings and export

### Creating Animations

**Animate time series (if per-timestep output is enabled):**

1. **File → Open** and select first VTK file
2. ParaView will detect file series (`.vt(0).vtu`, `.vt(1).vtu`, etc.)
3. **View → Animation View** (bottom panel appears)
4. In Animation View:
   - Configure start/end frames
   - Set frame rate (24 fps typical)
5. **Play** button to preview animation
6. **File → Export Animation** to save as video

### Modifying Color Maps

**Customize colors for better visualization:**

1. Right-click color bar in visualization
2. **Edit Color Map**
3. Choose from presets:
   - **Viridis:** Perceptually uniform
   - **Cool to Warm:** Blue-Red diverging
   - **Rainbow:** Traditional spectrum
   - **Grayscale:** Publication quality
4. Or define custom map by clicking points on the gradient

### Measuring Properties

**Interactive measurement tool:**

1. **Filters → Temporal Interpolator** (for time series)
2. Use **Ruler** tool (top toolbar) to measure distances
3. **Measure on Surface** to get scalar values at specific points

### Statistics and Analysis

**Compute model statistics:**

1. **Filters → Statistical Attributes**
2. Generates min, max, mean, standard deviation
3. View in **View → Statistics View**

---

## Troubleshooting

### File Won't Open

**Error: "Cannot find VTK file"**
- Ensure full path is correct
- Check that `save_vtk_files = True` in model parameters
- Verify model ran successfully (check console output)
- Look for `.vtu` files in `model_output/*/vtk_files/`

**Error: "Unsupported file format"**
- Ensure file is `.vtu` (XML unstructured grid), not `.vtu.bak` or other extension
- Verify file is not corrupted: `file model_output/scenario/vtk_files/*.vtu`

### Visualization Issues

**Colors don't show correctly**
- Right-click color bar → **Edit Color Map**
- Check that data range is appropriate
- Ensure variable contains non-constant values

**Mesh looks distorted**
- This is normal for unstructured grids
- Adjust viewing angle using mouse
- Use **View → Reset Camera** to reset perspective

**Performance is slow with large files**
- File may be too large for available RAM
- Reduce mesh resolution in model parameters
- Use **Filters → Downsample** to reduce mesh density
- Consider only visualizing final timestep (set `save_vtk_files_all_steps = False`)

### File Format Questions

**What's the difference between Elements and FaceElements?**
- **Elements:** Cell-centered data (primary visualization)
- **FaceElements:** Face-centered data (less commonly used)
- Use **Elements** files for standard visualization

**Why are files so large?**
- GroMPy uses 18 variables per mesh
- Files already use float32 precision and compression
- For 50×100 cell mesh: typically 1-2 MB (compressed)
- For larger meshes (>100×200 cells): may be 10-50 MB
- Consider reducing mesh resolution if files exceed available storage

### Getting Help

**ParaView Documentation:**
- Official tutorials: https://www.paraview.org/tutorials/
- Built-in help: **Help → Documentation** (in ParaView)

**GroMPy-couple Issues:**
- Check README.md in main repository
- Review CHANGELOG.md for VTK-related updates
- Open GitHub issue: https://github.com/ElcoLuijendijk/GroMPy-couple/issues

---

## Example Workflows

### Workflow 1: Studying Saltwater Intrusion

1. Run coastal aquifer model with `save_vtk_files = True`
2. Open final output VTK file in ParaView
3. **Visualize concentration:**
   - Color by "concentration"
   - Use Cool-Warm colormap (blue=fresh, red=salt)
4. **Show seawater region:**
   - Apply Threshold filter: `sea_surface_mask > 0.5`
5. **Add velocity field:**
   - Apply Glyph filter with arrows
   - Scale by flux_magnitude
6. **Export image:**
   - File → Export Scene as PNG

### Workflow 2: Analyzing Flow Paths in Terrestrial Aquifer

1. Open VTK file
2. **Visualize hydraulic head:**
   - Color by "hydraulic_head"
3. **Show flow direction:**
   - Apply Streamline filter
   - Seed from top surface (recharge area)
4. **Highlight recharge zone:**
   - Apply Threshold: `recharge_bnd > 0.5`
   - Make this layer semi-transparent
5. **Export for publication:**
   - Set high resolution (3000×2000 px)
   - Save as PDF for vector graphics

### Workflow 3: Examining Pressure-Concentration Coupling

1. Open VTK file
2. **Create split view:**
   - View → Split Horizontal
   - Load same VTK file in second view
3. **First view:** Color by pressure
4. **Second view:** Color by concentration
5. **Rotate both views identically** to see correlation
6. Useful for understanding coupling effects

---

## Tips and Best Practices

### General Tips

- **Start simple:** Open file, color by single variable first
- **Use meaningful colormaps:** Blue-red for diverging data, single color for sequential
- **Adjust opacity:** Use transparency to see overlapping features
- **Reset view:** Use **View → Reset Camera** if lost orientation
- **Use Undo:** Edit → Undo (Ctrl+Z) to revert unwanted changes

### For Publishing

- Export as **PDF** for presentations (vector graphics, clean text)
- Export as **PNG** with high resolution (3000×2000 px) for papers
- Use **grayscale colormap** for black & white printing
- Include **color bar** with units for clarity

### For Large Models

- Disable `save_vtk_files_all_steps` to reduce files
- Use smaller mesh resolution in model parameters
- Apply Threshold to focus on region of interest
- Use Slice for 2D visualization of 3D results

---

## Common Colormaps for Hydrogeology

| Variable | Recommended Colormap | Reasoning |
|----------|----------------------|-----------|
| Pressure | Cool to Warm | Shows low (cool) to high (warm) pressure |
| Concentration | Viridis | Perceptually uniform, good for printing |
| Velocity | Plasma | Shows slow (dark) to fast (bright) flow |
| Permeability | Viridis | Highlights heterogeneity clearly |
| Masks | Gray | Simple 0/1 visualization |

---

## Summary

ParaView is a powerful tool for visualizing GroMPy-couple results. Start with simple visualizations (single variable coloring), then add complexity (glyphs, streamlines, thresholds) as needed. Refer to this guide and ParaView's built-in help for specific tasks.

For questions about GroMPy-couple VTK output format, see `readme.md` VTK section.

**Happy visualizing!**
