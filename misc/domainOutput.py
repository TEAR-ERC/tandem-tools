#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:02:52 2024

@author: bar
"""

#%%

import numpy as np
import pyvista as pv
import meshio
import matplotlib.pyplot as plt
#%%



class mesh:
    def __init__(self, file_path):
        """
        Initialize the GmshMesh object.

        Parameters:
            file_path (str): Path to the GMSH mesh file (typically a .msh file).
        """
        self.file_path = file_path
        self.mesh = None
        self.physical_curves = {}      # phys_id -> list of line‐segment arrays
        self.phys_id_to_name = {}      # phys_id -> group name
        self.read_mesh()

    def read_mesh(self):
        """
        Reads the GMSH mesh file and extracts the 1D physical curves (physical groups).
        If none are defined, creates default IDs 1, 3, and 5.
        """
        self.mesh = meshio.read(self.file_path)

        # 1) Build mapping from physical id -> name for any 1D groups
        if self.mesh.field_data:
            for name, (phys_id, dim) in self.mesh.field_data.items():
                if dim == 1:
                    self.phys_id_to_name[phys_id] = name

        # 2) Extract any actual physical-curve data
        if "gmsh:physical" in self.mesh.cell_data:
            for cell_block, phys_data in zip(self.mesh.cells, 
                                             self.mesh.cell_data["gmsh:physical"]):
                if cell_block.type == "line":
                    for seg, phys in zip(cell_block.data, phys_data):
                        self.physical_curves.setdefault(phys, []).append(seg)

        # 3) Fallback defaults if _no_ physical curves were found
        if not self.physical_curves:
            defaults = {
                1: "surface",
                3: "fault",
                5: "direchelt"
            }
            for phys_id, name in defaults.items():
                # ensure the name mapping exists
                self.phys_id_to_name.setdefault(phys_id, name)
                # ensure the group exists (empty list if nothing to plot)
                self.physical_curves.setdefault(phys_id, [])

        return self.mesh, self.physical_curves, self.phys_id_to_name

    def plot_physical_curves(self, ax=None, show_triangles=True):
        """
        Plots the physical curves, or the default 1/3/5 groups if none were in the file.
        Optionally overlays triangles in light grey.
        """
        if self.mesh is None:
            raise ValueError("Mesh not loaded. Call read_mesh() first.")
        if ax is None:
            fig, ax = plt.subplots()

        # plot triangles
        if show_triangles:
            for cb in self.mesh.cells:
                if cb.type == "triangle":
                    tris = self.mesh.points[cb.data]
                    for tri in tris:
                        pts = tri[:, :2]
                        pts = np.vstack([pts, pts[0]])
                        ax.plot(pts[:,0], pts[:,1], color="lightgrey", linewidth=0.1,alpha=0.7)

        cmap = plt.get_cmap("tab10")
        for i, (phys, segs) in enumerate(self.physical_curves.items()):
            label = self.phys_id_to_name.get(phys, f"Physical {phys}")
            color = cmap(i % 10)
            for j, seg in enumerate(segs):
                pts = self.mesh.points[seg, :2]
                if j == 0:
                    ax.plot(pts[:,0], pts[:,1], color=color, linewidth=2, label=label)
                else:
                    ax.plot(pts[:,0], pts[:,1], color=color, linewidth=2)

        #ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        #ax.legend()
        #plt.title("Physical Curves with Mesh Triangles")
        plt.show()


#%%
# Define the StrainStressAnalyzer class with the updated compute_principal_stresses method
class domainFile:
    """
    Analyzes strain and stress data from a VTU file, computes principal stresses,
    and visualizes the results.

    Attributes:
        filename (str): Path to the VTU file.
        lambda_ (float or np.ndarray): Lame parameter λ.
        mu (float or np.ndarray): Lame parameter μ.
        mesh (pyvista.UnstructuredGrid): The loaded mesh from the VTU file.
    """
    
    def __init__(self, filename, lambda_=30e9, mu=30e9):
        """
        Initializes the StrainStressAnalyzer with the VTU file and Lame parameters.

        Args:
            filename (str): Path to the VTU file.
            lambda_ (float or np.ndarray): Lame parameter λ.
            mu (float or np.ndarray): Lame parameter μ.
        """
        self.filename = filename
        self.lambda_ = lambda_
        self.mu = mu
        self.mesh = pv.read(filename)
        print(f"VTU file '{filename}' loaded successfully.")
        self._extract_strain_derivatives()
        self.compute_all()
    def _extract_strain_derivatives(self):
        """Extracts strain derivative fields from the mesh."""
        required_fields = ["du0_d0", "du1_d1", "du1_d0", "du0_d1"]
        for field in required_fields:
            if field not in self.mesh.point_data:
                raise ValueError(f"Field '{field}' not found in the VTU file.")
            setattr(self, field, self.mesh.point_data[field])
            print(f"Strain derivative '{field}' extracted.")
    
    def compute_strain(self):
        """Computes strain components ε_xx, ε_yy, and ε_xy and stores them in the mesh."""
        self.mesh["Strain_E_xx"] = self.du0_d0
        self.mesh["Strain_E_yy"] = self.du1_d1
        self.mesh["Strain_E_xy"] = 0.5 * (self.du1_d0 + self.du0_d1)
        print("Strain components ε_xx, ε_yy, and ε_xy computed and stored in the mesh.")
    
    def compute_stress(self):
        """Computes stress components σ_xx, σ_yy, and σ_xy using plane strain relations."""
        # Handle scalar or array Lame parameters
        if np.isscalar(self.lambda_) and np.isscalar(self.mu):
            lambda_array = np.full_like(self.mesh["Strain_E_xx"], self.lambda_)
            mu_array = np.full_like(self.mesh["Strain_E_xx"], self.mu)
            print("Using global Lame parameters for stress computation.")
        else:
            lambda_array = self.lambda_
            mu_array = self.mu
            print("Using per-node Lame parameters for stress computation.")
        
        # Compute stress components
        self.mesh["Stress_S_xx"] = lambda_array * (self.mesh["Strain_E_xx"] + self.mesh["Strain_E_yy"]) + 2 * mu_array * self.mesh["Strain_E_xx"]
        self.mesh["Stress_S_yy"] = lambda_array * (self.mesh["Strain_E_xx"] + self.mesh["Strain_E_yy"]) + 2 * mu_array * self.mesh["Strain_E_yy"]
        self.mesh["Stress_S_xy"] = 2 * mu_array * self.mesh["Strain_E_xy"]
        print("Stress components σ_xx, σ_yy, and σ_xy computed and stored in the mesh.")
    
    def compute_principal_stresses(self):
        """
        Computes principal stresses and their directions, ensuring that σ1
        is the most compressive principal stress. The principal angle is
        calculated using arctan2 and adjusted to fall within -90° to +90°.
        
        The results are stored in the mesh's point data as:
            - "Principal_Sigma1": Most compressive principal stress (σ1)
            - "Principal_Sigma2": Less compressive principal stress (σ2)
            - "Principal_Angles": Orientation angle (θ_p) in degrees
        """
        sigma_xx = self.mesh["Stress_S_xx"]
        sigma_yy = self.mesh["Stress_S_yy"]
        sigma_xy = self.mesh["Stress_S_xy"]

        # Calculate average stress
        sigma_avg = 0.5 * (sigma_xx + sigma_yy)

        # Calculate the radius for principal stress calculation
        sigma_diff = 0.5 * (sigma_xx - sigma_yy)
        radius = np.sqrt(sigma_diff**2 + sigma_xy**2)

        # Compute initial principal stresses
        principal_sigma1_initial = sigma_avg + radius
        principal_sigma2_initial = sigma_avg - radius

        # Assign sigma1 as the most compressive stress (most negative)
        principal_sigma1 = np.minimum(principal_sigma1_initial, principal_sigma2_initial)
        principal_sigma2 = np.maximum(principal_sigma1_initial, principal_sigma2_initial)

        # Compute principal angle using arctan2 for correct quadrant determination
        theta_p_rad_full = np.arctan2(2 * sigma_xy, sigma_xx - sigma_yy)
        theta_p_full_deg = np.degrees(theta_p_rad_full)

        # Adjust angle to fall within -90° to +90°
        theta_p_deg_adjusted = np.where(
            theta_p_full_deg > 90,
            theta_p_full_deg - 180,
            np.where(theta_p_full_deg <= -90, theta_p_full_deg + 180, theta_p_full_deg)
        )

        # Compute the final principal angle by halving the adjusted angle
        theta_p_deg = 0.5 * theta_p_deg_adjusted

        # Store the results in the mesh's point data
        self.mesh["Principal_Sigma1"] = principal_sigma1  # Most compressive principal stress
        self.mesh["Principal_Sigma2"] = principal_sigma2  # Less compressive principal stress
        self.mesh["Principal_Angles"] = theta_p_deg       # Principal angle in degrees

        print("Principal stresses σ1 (most compressive), σ2 (less compressive), and angles computed and stored in the mesh.")
    
    def plot_quantity(self, quantity_name):
        """
        Plots the specified quantity from the mesh.

        Args:
            quantity_name (str): The name of the quantity to plot.
        """
        if quantity_name not in self.mesh.point_data and quantity_name not in self.mesh.cell_data:
            raise ValueError(f"Quantity '{quantity_name}' not found in the mesh.")
        
        plotter = pv.Plotter()
        plotter.add_title(f"{quantity_name} Visualization")
        plotter.add_mesh(self.mesh, scalars=quantity_name, cmap="bwr", show_edges=True)
        plotter.add_scalar_bar(title=quantity_name)
        plotter.view_xy()
        plotter.show_bounds(grid='front', location='outer', all_edges=True)
        plotter.enable_eye_dome_lighting()  # Enhances depth and shows coordinates under the mouse
    
        plotter.show()
        print(f"Plotted '{quantity_name}'.")
    
    def plot_principal_directions(self, scale=0.1):
        """
        Plots arrows indicating the direction of principal stresses (σ1) on the mesh.

        Args:
            scale (float): Scaling factor for the arrow lengths.
        """
        if "Principal_Angles" not in self.mesh.point_data:
            raise ValueError("Principal angles not found. Compute principal stresses first.")
        
        angles = self.mesh.point_data["Principal_Angles"]
        points = self.mesh.points
        vectors = np.column_stack((
            np.cos(np.radians(angles)),
            np.sin(np.radians(angles)),
            np.zeros_like(angles)
        ))
        vectors *= scale  # Scale for visibility
        
        # Add vectors to the mesh's point data
        self.mesh["Principal_Vectors"] = vectors
        
        # Create glyphs using the vectors
        glyphs = self.mesh.glyph(orient="Principal_Vectors", scale=False, factor=1.0)
        
        plotter = pv.Plotter()
        plotter.add_title("Principal Stress Directions (σ1)")
        plotter.add_mesh(
            self.mesh, 
            scalars="Stress_S_xx", 
            cmap="viridis", 
            opacity=0.6, 
            show_edges=True, 
            name="Stress Field"
        )
        plotter.add_mesh(
            glyphs, 
            color="red", 
            name="Principal Directions (σ1)"
        )
        plotter.add_scalar_bar(title="σ_xx (Pa)")
        plotter.add_legend([["Stress Field", "viridis"], ["Principal Directions (σ1)", "red"]])
        plotter.show()
        print("Principal stress directions (σ1) plotted successfully.")
    
    def compute_all(self):
        """Executes all computation steps: strain, stress, and principal stresses."""
        self.compute_strain()
        self.compute_stress()
        self.compute_principal_stresses()
        self.compute_coulomb_stress_change()
        print("All computations completed and stored in the mesh.")
    
    def plot_all(self):
        """
        Plots all relevant quantities:
        - Strain components ε_xx, ε_yy, ε_xy
        - Stress components σ_xx, σ_yy, σ_xy
        - Principal stresses σ1, σ2
        - Principal stress directions
        """
        quantities = ["Strain_E_xx", "Strain_E_yy", "Strain_E_xy",
                      "Stress_S_xx", "Stress_S_yy", "Stress_S_xy",
                      "Principal_Sigma1", "Principal_Sigma2"]
        
        for qty in quantities:
            self.plot_quantity(qty)
        
        self.plot_principal_directions(scale=0.1)
        print("All plots generated successfully.")
        
        
    def compute_coulomb_stress_change(self, dip_angle=30, friction=0.6):
        """
        Computes the column (Coulomb) stress change on a fault plane with a specified dip angle 
        (in degrees) and friction coefficient.

        The stress components (σ_nn, σ_nt) on the fault are computed as:
            σ_nn = σₓₓ·cos²θ + σᵧᵧ·sin²θ + 2σₓᵧ·sinθ·cosθ
            σ_nt = (σᵧᵧ - σₓₓ)·sinθ·cosθ + σₓᵧ·(cos²θ - sin²θ)
            
        Then, the column stress change is:
            Column Stress Change = σ_nt - (friction)·σ_nn

        The computed column stress change is stored in the mesh as "Column_Stress_Change"
        and printed.
        
        Args:
            dip_angle (float): Fault dip angle in degrees (default: 30).
            friction (float): Friction coefficient (default: 0.6).
        """
        # Ensure stress fields are computed
        if "Stress_S_xx" not in self.mesh.point_data or "Stress_S_yy" not in self.mesh.point_data or "Stress_S_xy" not in self.mesh.point_data:
            raise ValueError("Stress fields not computed. Please run compute_stress() first.")
        
        s_xx = self.mesh["Stress_S_xx"]
        s_yy = self.mesh["Stress_S_yy"]
        s_xy = self.mesh["Stress_S_xy"]
        
        theta = np.radians(dip_angle)
        c = np.cos(theta)
        s = np.sin(theta)
        
        sigma_nn = s_xx * c**2 + s_yy * s**2 + 2 * s_xy * s * c
        sigma_nt = (s_yy - s_xx) * s * c + s_xy * (c**2 - s**2)
        
        column_stress_change = sigma_nt - friction * sigma_nn
        self.mesh["Column_Stress_Change"] = column_stress_change
        

    
    
    def ExportUniformMeshForQunaitity(self, x, y,quantity_name):
        
        x,y=np.meshgrid(x,y)
        z=np.zeros_like(x)
        shape=x.shape
        rectGrid=pv.StructuredGrid(x,y,z)
        tandem_mesh_in_recet=rectGrid.sample(self.mesh)
        quantity=tandem_mesh_in_recet[quantity_name].reshape([shape[1],shape[0]]).T
        
        return quantity






