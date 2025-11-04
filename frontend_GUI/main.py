import sys
import PyQt5.QtWidgets as pyqtW
import PyQt5.QtGui as pyqtG
import PyQt5.QtCore as pyqtC

class MainWindow(pyqtW.QMainWindow): 

    def __init__(self): 
        super().__init__() 

        self.initUI()


    def initUI(self): 
        pass



if __name__ == "__main__": 

    app = pyqtW.QApplication(sys.argv)

    window = MainWindow() 
    window.show()
    sys.exit(app.exec_())

