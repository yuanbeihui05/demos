import numpy as np
import pyvista as pv

xyz_poly_coefs = (
                  np.array([
                    [[.5, -20.5, 7.875, 0, 0], [-15.0079, 28.0992, -72.75, 3.5, 0], [.0416667, -31.5825, .125, 0, 0], [0, -6.5, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, -14.5833, .916667, 0, 0], [-11.125, -2.72857, -33.4333, 1, 0], [-10.4893, 18.1516, 1.33333, 0, 0], [0, -46.4016, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[.75, -22.375, .9, 0, 0], [-14.2036, 10.5786, -51.6833, .8, 0], [-7.91429, -1.25, 1, 0, 0], [0, -39.6571, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[.5, -5.49206, 2.5, 0, 0], [-20.5, 28.6508, -22.0952, 8, 0], [7.875, .25, 2.5, 0, 0], [0, 28, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, -3.45833, 8.25, 0, 0], [-14.5833, 33.7286, -40.6667, 28, 0], [.916667, .6, 1, 0, 0], [0, 2.66667, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[.75, -8.17143, 1.15, 0, 0], [-22.375, 35.9714, -24.9429, 1, 0], [.9, 3.33333, .3, 0, 0], [0, 2.8, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[.5, 15.0079, 1.83333, 0, 0], [-5.49206, -12.6667, 18.6825, 6, 0], [2.5, 1.11111, .111111, 0, 0], [0, 2, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, 11.125, 1.88571, 0, 0], [-3.45833, -10.3119, 10.7071, 5.14286, 0], [8.25, .75, 1, 0, 0], [0, 5, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[.75, 14.2036, 1.75, 0, 0], [-8.17143, -17.3286, 18.0738, 4, 0], [1.15, 4, 1, 0, 0], [0, 3.6, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[.5, 20.5, 9.5, 0, 0], [15.0079, -13.6825, 1, 2, 0], [1.83333, 4, .8, 0, 0], [0, 1.33333, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, 14.5833, 6.375, 0, 0], [11.125, -8.16667, 1.33333, 1.5, 0], [1.88571, 4, 1, 0, 0], [0, 2.4, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[.75, 22.375, 10.5, 0, 0], [14.2036, -7.57857, .5, 2, 0], [1.75, 3.33333, 3.5, 0, 0], [0, 3, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[.5, 5.49206, 2.33333, 0, 0], [20.5, -8.31746, .222222, 1.33333, 0], [9.5, 2, 3, 0, 0], [0, 36, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, 3.45833, 1.33333, 0, 0], [14.5833, -8.25, 1.33333, 1.33333, 0], [6.375, 2.33333, 2.25, 0, 0], [0, 24, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[.75, 8.17143, 1.08571, 0, 0], [22.375, -16.1714, 2, 1.14286, 0], [10.5, 2.25, .777778, 0, 0], [0, 40, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[.5, -15.0079, .0416667, 0, 0], [5.49206, 2, 1.4, 6.66667, 0], [2.33333, 1.42857, 2, 0, 0], [0, 8, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1, -11.125, -10.4893, 0, 0], [3.45833, 4, 1.55556, 4.44444, 0], [1.33333, .25, 1.6, 0, 0], [0, 4, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[.75, -14.2036, -7.91429, 0, 0], [8.17143, 2, 4.5, 8, 0], [1.08571, 10, 1.25, 0, 0], [0, 3.2, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  )
uv_ranges = [[[0, 1],[0, 1]]]*6

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
