import sys
import PyQt5.QtWidgets as pyqtW
import PyQt5.QtGui as pyqtG
import PyQt5.QtCore as pyqtC

class MainWindow(pyqtW.QMainWindow): 

    def __init__(self): 
        super().__init__() 

        self.initUI()


    def initUI(self): 
        
        self.setWindowTitle("MolecuSim")
        layout = pyqtW.QGridLayout()

        nav_bar = pyqtW.QMenuBar()
        nav_bar.setStyleSheet("background: lightgrey; width: 1000px; max-height: 50px;")
        side_bar = pyqtW.QMenuBar()
        side_bar.setStyleSheet("background: lightgrey; max-width: 100px; height: 1200px")
        layout.addWidget(nav_bar, 0, 1)
        layout.addWidget(side_bar, 0, 0, 2, 1)

        main_window = pyqtW.QMenuBar()
        main_window.setStyleSheet("background: white; width: 1000px; height: 1000px;")
        layout.addWidget(main_window, 1, 1)

        widget = pyqtW.QWidget()
        widget.setLayout(layout)

        self.setCentralWidget(widget)



if __name__ == "__main__": 

    app = pyqtW.QApplication(sys.argv)

    window = MainWindow() 
    window.show()
    sys.exit(app.exec_())

