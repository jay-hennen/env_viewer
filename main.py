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
import matplotlib


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
        
        self.names_view.itemClicked.connect(self.changefig)
        f = Figure()
        self.canvas = FigureCanvas(f)
        self.cur_fig = self.fig_dict['default'] = FigManager('default', self.canvas, f)
        self.mplvl.addWidget(self.canvas)

        self.toolbar = NavigationToolbar(self.canvas,
                                         self.mplwindow,
                                         coordinates=True)
        self.mplvl.addWidget(self.toolbar)
        
        self.init_sliders()

    def init_sliders(self):
        self.slider_data = {'cur_time':self.date_time_box.dateTime(),
                            'time_idx_max':self.time_step_slider.maximum,
                            'depth_idx_max':self.depth_slider.maximum}
#         self.button_open_folder.clicked.connect(self.open_netCDF)
        self.time_step_slider._t_idx = 0
        self.time_step_slider.valueChanged.connect(self._time_changed)
        self.time_step_slider.sliderReleased.connect(self._time_changed)

    def _time_changed(self, idx=None):
        if not self.time_step_slider.isSliderDown():
            if idx is None:
                idx = self.time_step_slider._t_idx
            self.changefig(None, idx)
        else:
            self.time_step_slider._t_idx = idx

    def setup_sliders(self, env_obj):
        print 'slider max = {0}'.format(len(env_obj.time.time) - 1)
        self.time_step_slider.setMaximum(len(env_obj.time.time) - 1)
        self.time_step_slider.setValue(0)

    def open_netCDF(self):
        self.names_view.clear()
        self.fig_dict.clear()
        filename = QtWidgets.QFileDialog.getOpenFileName(self, "Pick a .nc")[0]

        if filename is not None:
            from gnome.environment import env_from_netCDF
            try:
                env = env_from_netCDF(filename)
             
                self.names_view_data = {}
                for e in env:
                    f = Figure()
                    mpl_objs = self.plot_env_obj(e, f)

                    self.addfig(e.name, f, e, mpl_objs)
            except Exception:
                raise

    def plot_env_obj(self, e, fig=None, time=None, time_idx=None, mpl_objs={}):
        if fig is None:
            fig = self.cur_fig.fig
        plate = ccrs.PlateCarree()
        pole = projection = ccrs.NorthPolarStereo()
        lon = e.grid.node_lon
        lat = e.grid.node_lat
        scl = (lon.shape[0] / 100 + 1, lon.shape[1] / 100 + 1)
        print scl
        lon = lon[::scl[0], ::scl[1]]
        lat = lat[::scl[0], ::scl[1]]
        pts = np.column_stack((lon.reshape(-1), lat.reshape(-1), np.zeros_like(lat.reshape(-1))))
        lon = lon - 360
        if time is None:
            if time_idx is None:
                t = e.time.min_time
            else:
                t = e.time.time[time_idx]
        fig_children = fig.get_children()
        if len(fig_children) > 1:
             p = fig_children[1]
        else:
            p = fig.add_subplot(111, projection=pole)
            mpl_objs = {'plot':p}
            p.coastlines('50m')
            p.set_extent([-180, 180, lat.min(), 90], ccrs.PlateCarree())
            p.add_feature(cartopy.feature.OCEAN, zorder=0)
            p.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
            p.gridlines()
#                         r = p.quiver(lon, lat, u, v, scale=None, units='xy', scale_units='xy', width=0.030)
        if isinstance(e, (GridVectorProp,)):
            vels = e.at(pts, time=t, interpolation='linear')
            u = vels[:, 0]
            v = vels[:, 1]
            u = np.ma.masked_equal(u, 0.0)
            v = np.ma.masked_equal(v, 0.0)
            u, v = self.convert_uv_to_delta(lat.reshape(-1), (u, v))
            u = u.reshape(lon.shape)
            v = v.reshape(lat.shape)
            if (time is not None or time_idx is not None) and 'quiver' in mpl_objs:
                q = mpl_objs['quiver']
                q.set_UVC(u, v)
                p.draw_artist(q)
            else:
                q = p.quiver(lon, lat, u, v, transform=ccrs.PlateCarree(), picker=5)
                mpl_objs['quiver'] = q
        if isinstance(e, GriddedProp):
            vals = e.at(pts, t, interpolation='linear')
            cb = None
            if (time is not None or time_idx is not None):
                cb = fig.axes[1]
                cb.clear()
            c = p.contourf(lon, lat, vals.reshape(lon.shape), transform=plate)
            b = fig.colorbar(c, cax=cb)
#                 mpl_objs['contour'] = c
#                 p._c = c
        return mpl_objs

    def changefig(self, item=None, time_idx=None):
        fm = None
        if item is not None:
            # Change of item, so resize/redraw everything
            text = item.text()
            fm = self.fig_dict[text]
            self.cur_fig = fm
            env_obj = fm.env_obj
            if env_obj is not None:
                self.setup_sliders(env_obj)
    
            fm.fig.set_size_inches(fm.width / fm.fig.dpi, fm.height / fm.fig.dpi)
        if time_idx is not None:
            if fm is None:
                fm = self.cur_fig
            if fm.env_obj is not None:
                self.plot_env_obj(fm.env_obj, fm.fig, time_idx=time_idx, mpl_objs=fm.mpl_objs)
            print 'Time index to {0}'.format(time_idx)
            # change displayed data on fm
            
            
        self.canvas.figure = fm.fig
        self.canvas.draw()
            # Time change, so redraw current item with new data

    def addfig(self, name, fig, env_obj=None, mpl_objs=None):
        self.fig_dict[name] = FigManager(name, self.canvas, fig, env_obj, mpl_objs)
        self.names_view.addItem(copy.deepcopy(name))

    def convert_uv_to_delta(self, lat, (u, v)):
        scale = 100
        timestep = 3600
        u *= scale * 8.9992801e-06
        v *= scale * 8.9992801e-06
        u /= np.cos(np.deg2rad(lat))
        return (u, v)

#     def change_bb(self, ax):

class FigManager(object):
    def __init__(self,
                 name,
                 canvas,
                 fig=None,
                 env_obj=None,
                 mpl_objs=None,
                 datetimeedit=None,
                 time_step_slider=None
                 ):
        self.name = name
        if fig is None:
            fig = Figure()
        self.fig = fig
        self.fig.canvas = canvas
        self.canvas = canvas
        # below are state variables for what is being displayed
        self.env_obj = env_obj
        self.mpl_objs = mpl_objs
        self.datetimeedit = datetimeedit
        self.time_step_slider = None
    
    @property
    def width(self):
        return self.canvas.size().width()
    
    @property
    def height(self):
        return self.canvas.size().height()
    
    @property
    def datetime_index(self):
        if self.env_obj is None:
            raise NotImplementedError("not available for non-env-objects")
        else:
            return e.time.index_of(self.datetimeedit.datetime().toPyDateTime())

    @property
    def time_step_index(self):
        if self.env_obj is None:
            raise NotImplementedError("not available for non-env-objects")
        else:
            return self.time_step_slider.value

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
