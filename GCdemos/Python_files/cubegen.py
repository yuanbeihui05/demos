import numpy as np
import pyvista as pv

xyz_poly_coefs = (
                  np.array([
                    [[417601/3, -43200, -45760, 0, 0], [-19200, 0, 3840, 0, 0], [-53200, 8640, 31440, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[622081/8, -12960, -23040, 0, 0], [-100800, 0, 20160, 0, 0], [-35928, 2592, 18792, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[149404, -16200, -45960, 0, 0], [-14400, 0, 2880, 0, 0], [-427020/7, 3240, 237960/7, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[-57599/3, 46400, 39600, 0, 0], [-43200, 0, 8640, 0, 0], [17040, -9280, -18240, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[-1266047/8, 25920, 73368, 0, 0], [-12960, 0, 2592, 0, 0], [63288, -5184, -42192, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[-133712/7, 49800, 303120/7, 0, 0], [-16200, 0, 3240, 0, 0], [155700/7, -9960, -149580/7, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[-120959/3, 46400, 44160, 0, 0], [-19200, 0, 3840, 0, 0], [14160, -9280, -18000, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[-5759/8, 25920, 16560, 0, 0], [-100800, 0, 20160, 0, 0], [-6408, -5184, -2952, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[8284, 49800, 29160, 0, 0], [-14400, 0, 2880, 0, 0], [-48600/7, -9960, -48060/7, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[115201/3, 46400, 20400, 0, 0], [-43200, 0, 8640, 0, 0], [-2160, -9280, -6720, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1153153/8, 25920, -27432, 0, 0], [-12960, 0, 2592, 0, 0], [-37512, -5184, 18288, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[168688/7, 49800, 202320/7, 0, 0], [-16200, 0, 3240, 0, 0], [54900/7, -9960, -89100/7, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[267841/3, 46400, 960, 0, 0], [-19200, 0, 3840, 0, 0], [-29040, -9280, 7920, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[305281/8, 25920, 3600, 0, 0], [-100800, 0, 20160, 0, 0], [-19368, -5184, 4824, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[56884, 49800, 12960, 0, 0], [-14400, 0, 2880, 0, 0], [-162000/7, -9960, 19980/7, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    ]).transpose(1,2,0),

                  np.array([
                    [[1/3, -43200, 640, 0, 0], [-19200, 0, 3840, 0, 0], [-6800, 8640, 3600, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[1/8, -12960, 2880, 0, 0], [-100800, 0, 20160, 0, 0], [-10008, 2592, 3240, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
                    [[4, -16200, 3840, 0, 0], [-14400, 0, 2880, 0, 0], [-78420/7, 3240, 28800/7, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
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
