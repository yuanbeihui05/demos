import numpy as np
import pyvista as pv

xyz_poly_coefs = (
                  np.array([
                    [[5/6, -267/35, -97/70, 0, 0], [232/35, -472/35, 4/5, 18/7, 0], [5, 5, 3/4, 0, 0], [9/10, 1/3, 0, 0, 0], [3/2, 0, 0, 0, 0]],
                    [[1/4, -44/3, -10/3, 0, 0], [43/3, -571/24, 1/3, 8/3, 0], [5/3, 1/2, 4/5, 0, 0], [10/7, 10/3, 0, 0, 0], [5/2, 0, 0, 0, 0]],
                    [[5/8, -1439/45, -67/9, 0, 0], [1412/45, -2213/45, 1/5, 14, 0], [1, 5/3, 6, 0, 0], [2, 2/3, 0, 0, 0], [9/2, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[5/6, 1, 3/4, 2/3, 1/3], [-267/35, 1, 1, 3/10, 0], [-97/70, 9/7, 1/6, 0, 0], [0, 1/5, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1/4, 1/3, 8/9, 9/2, 10/3], [-44/3, 9/8, 3/2, 8/3, 0], [-10/3, 1, 5/3, 0, 0], [0, 4, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[5/8, 3/5, 7/6, 7/5, 8/5], [-1439/45, 1/9, 3/5, 1/2, 0], [-67/9, 2, 1/6, 0, 0], [0, 8/9, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  )
uv_ranges = [[[0, 1],[0, 1]]]*2

if __name__ == '__main__':
    pv.set_plot_theme('paraview')
    plotter = pv.Plotter()
    plotter.set_color_cycler([
        "#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", 
        "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"
        ])
    plotter.show_axes_all()
    plotter.show_grid()
    plotter.enable_anti_aliasing()

    for (u_range,v_range),coefs in zip(uv_ranges, xyz_poly_coefs): 
        u = np.linspace(u_range[0], u_range[1], 20)
        v = np.linspace(v_range[0], v_range[1], 20)
        value = np.polynomial.polynomial.polygrid2d(u, v, coefs)

        mesh = pv.StructuredGrid(value[0], value[1], value[2])   
        plotter.add_mesh(mesh, smooth_shading=True, split_sharp_edges=True)

    plotter.show()
