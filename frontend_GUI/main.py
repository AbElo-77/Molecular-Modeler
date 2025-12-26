import sys
import PyQt5.QtWidgets as pyqtW
import PyQt5.QtGui as pyqtG
import PyQt5.QtCore as pyqtC
import pyqtgraph as pg
import pyqtgraph.opengl as og

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

class MainWindow(pyqtW.QMainWindow): 

    def __init__(self): 
        super().__init__() 

        self.initUI()


    def initUI(self): 
        
        self.setWindowTitle("MolecuSim")

        central = pyqtW.QWidget()
        layout = pyqtW.QGridLayout()
        central.setLayout(layout)
        self.setCentralWidget(central)

        nav_bar = pyqtW.QFrame()
        nav_bar.setStyleSheet("background: lightgrey; max-height: 50px;")
        layout.addWidget(nav_bar, 0, 1)

        side_bar = pyqtW.QFrame()
        side_bar.setStyleSheet("background: lightgrey; max-width: 100px;")
        layout.addWidget(side_bar, 0, 0, 2, 1)

        gl_widget = og.GLViewWidget()
        gl_widget.setBackgroundColor('w')
        gl_widget.opts['distance'] = 20

        size = 500
        spacing = 50

        grid_xy = og.GLGridItem()
        grid_xy.setSize(size, size, 1)
        grid_xy.setColor((1, 1, 1, 1))
        grid_xy.setSpacing(spacing, spacing, 1)
        grid_xy.rotate(90, 1, 0, 0)
        grid_xy.setGLOptions('additive')
        gl_widget.addItem(grid_xy)

        grid_xz = og.GLGridItem()
        grid_xz.setSize(size, size, 1)
        grid_xz.setColor((0, 0, 0, 1))
        grid_xz.setSpacing(spacing, spacing, 1)
        grid_xz.setGLOptions('additive')
        gl_widget.addItem(grid_xz)

        grid_yz = og.GLGridItem()
        grid_yz.setSize(size, size, 1)
        grid_yz.setColor((0, 0, 0, 1))
        grid_yz.setSpacing(spacing, spacing, 1)
        grid_yz.setGLOptions('additive')
        grid_yz.rotate(90, 0, 1, 0)
        gl_widget.addItem(grid_yz)

        axis = og.GLAxisItem()
        axis.setGLOptions('opaque')
        sphere = og.MeshData.sphere(20, 40, 3)
        sphere_item = og.GLMeshItem(meshdata=sphere, smooth=True, color=(0, 0, 0, 0), shader='balloon')
        sphere_item.translate(1, 2, 3)
        
        gl_widget.addItem(axis)
        gl_widget.addItem(sphere_item)

        layout.addWidget(gl_widget, 1, 1)



if __name__ == "__main__": 

    app = pyqtW.QApplication(sys.argv)

    window = MainWindow() 
    window.show()
    sys.exit(app.exec_())

