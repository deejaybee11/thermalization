# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:40:08 2016

@author: dbro184
"""

import sys
from PyQt4 import QtCore, QtGui
from PyQt4.uic import loadUiType
import astropy.io.fits as fits
from matplotlib.figure import Figure
import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

Ui_MainWindow, QMainWindow = loadUiType('test.ui')

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        self.matplotLayout.addWidget(self.canvas)
        self.matplotLayout.addWidget(self.mpl_toolbar)
        self.mydata = fits.open("../../C++/thermalization/testexpand3D.fits")
        self.mydata = self.mydata[0].data
        self.sbZ.setMaximum(len(self.mydata[:,:,0]))
        self.sbZ.valueChanged.connect(self.show_img_spin)   
        self.show_img_spin() 

    def show_img_spin(self, ):

        self.ax.imshow(self.mydata[:,:,self.sbZ.value()], vmin=0, vmax = np.amax(self.mydata))        
        self.canvas.draw()
        self.ax.clear()

        
if __name__ == "__main__":
    
    app = QtGui.QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    myapp = Main()
    myapp.show()
    app.exec_()
#    sys.exit(0)
