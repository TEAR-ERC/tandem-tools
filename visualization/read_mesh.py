import meshio
import matplotlib.pyplot as plt
import numpy as np
'''
Originally written by Bar Oryan
Modified by Jeena Yun
'''

DEFAULT_GROUPS = {1: "surface", 3: "fault", 5: "dirichlet"}

class GmshMesh2D:
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
        self.dim=2
        self.read_mesh()

    def read_mesh(self) -> None:
        """Read the file and fill .mesh / .physical_curves / .phys_id_to_name"""
        self.mesh = meshio.read(self.file_path)

        # 1) map phys-ID → name  (only for 1-D groups)
        if self.mesh.field_data:
            for name, (phys_id, dim) in self.mesh.field_data.items():
                if dim == 1:                 # 1-D = curve
                    self.phys_id_to_name[phys_id] = name

        # 2) extract line elements that carry a physical tag
        if "gmsh:physical" in self.mesh.cell_data:
            for cb, phys in zip(self.mesh.cells,
                                self.mesh.cell_data["gmsh:physical"]):
                if cb.type == "line":
                    for seg, tag in zip(cb.data, phys):
                        self.physical_curves.setdefault(tag, []).append(seg)

        # 3) ensure default IDs exist (possibly empty lists)
        if not self.physical_curves:           # nothing found → fallback groups
            for pid, name in DEFAULT_GROUPS.items():
                self.phys_id_to_name.setdefault(pid, name)
                self.physical_curves.setdefault(pid, [])

    def plot_physical_curves(self, ax=None, show_triangles=True, save_dir=None):
        """
        Plots the physical curves, or the default 1/3/5 groups if none were in the file.
        Optionally overlays triangles in light grey.
        """
        if self.mesh is None:
            raise ValueError("Mesh not loaded. Call read_mesh() first.")
        if ax is None:
            import myplots
            mp = myplots.Figpref()
            plt.rcParams['font.size'] = '11'
            fig,ax = plt.subplots(figsize=(10,5))

        # plot triangles
        xmin, xmax = 1e10, -1e10
        ymin, ymax = 1e10, -1e10
        if show_triangles:
            for cb in self.mesh.cells:
                if cb.type == "triangle":
                    tris = self.mesh.points[cb.data]
                    for tri in tris:
                        pts = tri[:, :2]
                        pts = np.vstack([pts, pts[0]])
                        ax.plot(pts[:,0], pts[:,1], color="lightgrey", linewidth=0.25,alpha=0.7)
                        xmin, xmax = self.update_min_max(xmin, xmax, pts[:,0])
                        ymin, ymax = self.update_min_max(ymin, ymax, pts[:,1])

        # cmap = plt.get_cmap("tab10")
        for i, (phys, segs) in enumerate(self.physical_curves.items()):
            # label = self.phys_id_to_name.get(phys, f"Physical {phys}")
            # color = cmap(i % 10)
            if phys == 1:
                color = mp.myblue
                label = 'Free surface'
            elif phys == 3:
                color = 'k'
                label = 'Fault (rate-and-state friction)'
            elif phys == 5:
                color = mp.mypink
                label = 'Dirichlet'
            for j, seg in enumerate(segs):
                pts = self.mesh.points[seg, :2]
                if j == 0:
                    ax.plot(pts[:,0], pts[:,1], color=color, linewidth=2, label=label)
                else:
                    ax.plot(pts[:,0], pts[:,1], color=color, linewidth=2)
                xmin, xmax = self.update_min_max(xmin, xmax, pts[:,0])
                ymin, ymax = self.update_min_max(ymin, ymax, pts[:,1])

        #ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel("X [km]", fontsize=13)
        ax.set_ylabel("Y [km]", fontsize=13)
        xl = xmax - xmin
        yl = ymax - ymin
        ax.set_xlim(xmin-xl*0.01, xmax+xl*0.01)
        ax.set_ylim(ymin-yl*0.01, ymax+yl*0.01)
        ax.legend(fontsize=11, loc='lower right', bbox_to_anchor=(0.98, 0.02), bbox_transform=ax.transAxes)
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        #plt.title("Physical Curves with Mesh Triangles")
        plt.tight_layout()
        if save_dir is not None:
            plt.savefig('%s/mesh_BC.png'%(save_dir),dpi=300)
        
    def get_curve_points(self, phys_id, dim=None):
        """
        Extracts the endpoint coordinates of every line‐segment in a physical curve.
    
        Parameters
        ----------
        mesh : meshio.Mesh
            A meshio mesh object with mesh.points of shape (N,3).
        physical_curves : dict[int, list[array]]
            Mapping each physical‐curve ID to a list of 1D int‐arrays (node indices).
        phys_id : int
            The physical‐curve ID whose segments you want.
        dim : {2,3}, default 2
            Number of coordinates to return per node (2 → XY only, 3 → XYZ).
    
        Returns
        -------
        List[np.ndarray]
            A list of length-2 arrays of shape (2, dim), one per segment.
        """
        if dim is None:
            dim=self.dim
        #pts_list = []
        # grab the list of segments (each seg is e.g. array([n0, n1], dtype=int))
        pts_list: list[np.ndarray] = []
    
        for element in self.physical_curves.get(phys_id, []):
            # mesh.points[element] has shape (Nnodes, 3); slice → (Nnodes, dim)
            pts_list.append(self.mesh.points[element, :dim])
    
        return pts_list
    
    def update_min_max(self, current_min, current_max, val):
        new_min = current_min
        new_max = current_max
        if current_min > np.min(val):
            new_min = np.min(val)
        if current_max < np.max(val):
            new_max = np.max(val)
        return new_min, new_max