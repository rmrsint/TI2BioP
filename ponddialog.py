__author__ = 'reymolina'

import sys
from PyQt5 import QtGui, QtCore, QtWidgets

from ponderation import *
import shelve



class DialogPond(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.ui = Ui_Pondx()
        self.ui.setupUi(self)
        self.ui.Select_all.clicked.connect(self.SelectAll)
        self.ui.Clear_all.clicked.connect(self.ClearAll)
        self.ui.buttonBox.accepted.connect(self.accept)

    def Addponderations(self):
        names = shelve.open('datapond.db', flag='r')

        for i in range(4, len(names['Dname'])):
            item = QtWidgets.QListWidgetItem()
            item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable|QtCore.Qt.ItemIsEnabled)
            item.setCheckState(QtCore.Qt.Unchecked)
            item.setText(names['Dname'][i])
            self.ui.listWidget.addItem(item)
        names.close()

    def SelectAll(self):
        for i in range(self.ui.listWidget.count()):
            item = self.ui.listWidget.item(i)
            item.setCheckState(QtCore.Qt.Checked)

    def ClearAll(self):
        for i in range(self.ui.listWidget.count()):
            item = self.ui.listWidget.item(i)
            item.setCheckState(QtCore.Qt.Unchecked)

    def getvalues(self):
        selected = []
        for i in range(self.ui.listWidget.count()):
            item = QtWidgets.QListWidgetItem( self.ui.listWidget.item(i))
            if item.checkState() == QtCore.Qt.Checked:
                selected.append(i)
        return selected



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    myapp = DialogPond()
    myapp.Addponderations()
    if myapp.exec_():
        values = myapp.getvalues()
        print (values)


    #myapp.show()

    #sys.exit(app.exec_())








