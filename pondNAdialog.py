__author__ = 'reymolina'

import sys
from PyQt5 import QtGui, QtCore, QtWidgets

from ponderation_NA import *


class DialogPondRNA(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self,parent)
        self.ui=Ui_PondRNA()
        self.ui.setupUi(self)
        self.ui.Select_all.clicked.connect(self.SelectAll)
        self.ui.Clear_all.clicked.connect(self.ClearAll)
        self.ui.buttonBox.accepted.connect(self.accept)


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
            item = QtWidgets.QListWidgetItem(self.ui.listWidget.item(i))
            if item.checkState() == QtCore.Qt.Checked:
                selected.append(i)
        return selected


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    myapp = DialogPondRNA()
    if myapp.exec_():
        values = myapp.getvalues()
        print(values)
