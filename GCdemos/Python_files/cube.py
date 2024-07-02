import numpy as np
import pyvista as pv

xyz_poly_coefs = (
                  np.array([
                    [[0, 0, 0, 0, 0], [2/3, 0, -2/15, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 2/3, 0, 0, 0], [0, 0, 0, 0, 0], [0, -2/15, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, 0, -1/3, 0, 0], [0, 0, 0, 0, 0], [-1/3, 0, 1/5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[1, 0, -1/3, 0, 0], [0, 0, 0, 0, 0], [-1/3, 0, 1/5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 0, 0, 0, 0], [2/3, 0, -2/15, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 2/3, 0, 0, 0], [0, 0, 0, 0, 0], [0, -2/15, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[0, 0, 0, 0, 0], [2/3, 0, -2/15, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, 0, -1/3, 0, 0], [0, 0, 0, 0, 0], [-1/3, 0, 1/5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 2/3, 0, 0, 0], [0, 0, 0, 0, 0], [0, -2/15, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[-1, 0, 1/3, 0, 0], [0, 0, 0, 0, 0], [1/3, 0, -1/5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 0, 0, 0, 0], [2/3, 0, -2/15, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 2/3, 0, 0, 0], [0, 0, 0, 0, 0], [0, -2/15, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[0, 0, 0, 0, 0], [2/3, 0, -2/15, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[-1, 0, 1/3, 0, 0], [0, 0, 0, 0, 0], [1/3, 0, -1/5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 2/3, 0, 0, 0], [0, 0, 0, 0, 0], [0, -2/15, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[0, 0, 0, 0, 0], [2/3, 0, -2/15, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[0, 2/3, 0, 0, 0], [0, 0, 0, 0, 0], [0, -2/15, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[-1, 0, 1/3, 0, 0], [0, 0, 0, 0, 0], [1/3, 0, -1/5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  )
uv_ranges = [[[-1, 1],[-1, 1]]]*6

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
