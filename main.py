from PyQt5 import QtCore, QtGui, QtWidgets

import sys, os  # We need sys so that we can pass argv to QApplication
import copy
import env_viewer  # This file holds our MainWindow and all design related things
              # it also keeps events etc that we defined in Qt Designer

from gnome.environment.grid_property import GridVectorProp, GriddedProp

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from gnome.environment.grid_property import GridVectorProp

import cartopy
import cartopy.crs as ccrs


class EnvViewer(QtWidgets.QMainWindow, env_viewer.Ui_MainWindow):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
        self.button_open_folder.clicked.connect(self.open_netCDF)
                            # It sets up layout and widgets that are defined
        self.fig_dict = {}
        self.names_view_data = {}
        
        self.names_view.itemClicked.connect(self.changefig)
        
        fig = Figure()
        self.addmpl(fig)
        
    def open_netCDF(self):
        self.names_view.clear()
        self.names_view_data.clear()
        self.fig_dict.clear()
        self.env_objs = None
        filename = QtWidgets.QFileDialog.getOpenFileName(self, "Pick a .nc")[0]
 
        if filename is not None:
            from gnome.environment import env_from_netCDF
            try:
                env = env_from_netCDF(filename)
                self.env_objs = env
             
                self.names_view_data = {}
                for e in env:
                    self.names_view.addItem(e.name)
                    f = Figure()
                    plate = ccrs.PlateCarree()
                    pole = projection = ccrs.NorthPolarStereo()
                    self.names_view_data[e.name] = e
                    lon = e.grid.node_lon
                    lat = e.grid.node_lat
                    scl = (lon.shape[0] / 100 + 1, lon.shape[1] / 100 + 1)
                    print scl
                    lon = lon[::scl[0], ::scl[1]]
                    lat = lat[::scl[0], ::scl[1]]
                    pts = np.column_stack((lon.reshape(-1), lat.reshape(-1), np.zeros_like(lat.reshape(-1))))
                    lon = lon - 360
                    t = e.time.min_time
                    p = f.add_subplot(111, projection=pole)
                    p.coastlines('50m')
                    p.set_extent([-180, 180, lat.min(), 90], ccrs.PlateCarree())
                    p.add_feature(cartopy.feature.OCEAN, zorder=0)
                    p.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
                    p.gridlines()
#                         r = p.quiver(lon, lat, u, v, scale=None, units='xy', scale_units='xy', width=0.030)
                    if isinstance(e, (GridVectorProp,)):
                        vels = e.at(pts, t, interpolation='linear')
                        u = vels[:, 0]
                        v = vels[:, 1]
                        u = np.ma.masked_equal(u, 0.0)
                        v = np.ma.masked_equal(v, 0.0)
                        u, v = self.convert_uv_to_delta(lat.reshape(-1), (u, v))
                        u = u.reshape(lon.shape)
                        v = v.reshape(lat.shape)
                        r = p.quiver(lon, lat, u, v, transform=ccrs.PlateCarree())
                    if isinstance(e, GriddedProp):
                        vals = e.at(pts, t, interpolation='linear')
                        c = p.contourf(lon, lat, vals.reshape(lon.shape), transform=plate)
                        f.colorbar(c)
                    self.fig_dict[e.name] = f
            except Exception:
                raise

    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
                                         self.mplwindow,
                                         coordinates=True)
        self.mplvl.addWidget(self.toolbar)
#         self.addToolBar(self.toolbar)

    def rmmpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()

    def changefig(self, item):
        text = item.text()
        self.rmmpl()
        self.addmpl(self.fig_dict[text])

    def addfig(self, name, fig):
        self.fig_dict[name] = fig
        self.names_view.addItem(copy.deepcopy(name))

    def convert_uv_to_delta(self, lat, (u, v)):
        scale = 100
        timestep = 3600
        u *= scale * 8.9992801e-06
        v *= scale * 8.9992801e-06
        u /= np.cos(np.deg2rad(lat))
        return (u, v)

#     def change_bb(self, ax):

if __name__ == '__main__':  # if we're running file directly and not importing it
    import sys
    import numpy as np
        
    fig1 = Figure()
    ax1f1 = fig1.add_subplot(111)
    ax1f1.plot(np.random.rand(5))
 
    fig2 = Figure()
    ax1f2 = fig2.add_subplot(121)
    ax1f2.plot(np.random.rand(5))
    ax2f2 = fig2.add_subplot(122)
    ax2f2.plot(np.random.rand(10))
 
    fig3 = Figure()
    ax1f3 = fig3.add_subplot(111)
    ax1f3.pcolormesh(np.random.rand(20, 20))
 
    app = QtWidgets.QApplication(sys.argv)  # A new instance of QApplication
    main = EnvViewer()  # We set the form to be our ExampleApp (design)
    main.addfig('One plot', fig1)
    main.addfig('Two plots', fig2)
    main.addfig('Pcolormesh', fig3)
    main.show()
    app.exec_()
