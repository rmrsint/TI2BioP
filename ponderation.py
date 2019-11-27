# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ponderation.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Pondx(object):
    def setupUi(self, Pondx):
        Pondx.setObjectName("Pondx")
        Pondx.resize(568, 279)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/newPrefix/RNA_2_Icon.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Pondx.setWindowIcon(icon)
        self.verticalLayout = QtWidgets.QVBoxLayout(Pondx)
        self.verticalLayout.setObjectName("verticalLayout")
        self.listWidget = QtWidgets.QListWidget(Pondx)
        self.listWidget.setObjectName("listWidget")
        item = QtWidgets.QListWidgetItem()
        item.setCheckState(QtCore.Qt.Checked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        item.setFlags(QtCore.Qt.ItemIsUserCheckable|QtCore.Qt.ItemIsEnabled)
        item.setCheckState(QtCore.Qt.Unchecked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable|QtCore.Qt.ItemIsEnabled)
        item.setCheckState(QtCore.Qt.Unchecked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.listWidget.addItem(item)
        self.verticalLayout.addWidget(self.listWidget)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.Select_all = QtWidgets.QPushButton(Pondx)
        self.Select_all.setObjectName("Select_all")
        self.horizontalLayout.addWidget(self.Select_all)
        self.Clear_all = QtWidgets.QPushButton(Pondx)
        self.Clear_all.setObjectName("Clear_all")
        self.horizontalLayout.addWidget(self.Clear_all)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        self.buttonBox = QtWidgets.QDialogButtonBox(Pondx)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.horizontalLayout_2.addWidget(self.buttonBox)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(Pondx)
        self.buttonBox.accepted.connect(Pondx.accept)
        self.buttonBox.rejected.connect(Pondx.reject)
        QtCore.QMetaObject.connectSlotsByName(Pondx)

    def retranslateUi(self, Pondx):
        _translate = QtCore.QCoreApplication.translate
        Pondx.setWindowTitle(_translate("Pondx", "Ponderations"))
        __sortingEnabled = self.listWidget.isSortingEnabled()
        self.listWidget.setSortingEnabled(False)
        item = self.listWidget.item(0)
        item.setText(_translate("Pondx", "sum-Amber95"))
        item = self.listWidget.item(1)
        item.setText(_translate("Pondx", "sum-FCPAC free"))
        item = self.listWidget.item(2)
        item.setText(_translate("Pondx", "sum-FCPAC buried"))
        item = self.listWidget.item(3)
        item.setText(_translate("Pondx", "ISA"))
        self.listWidget.setSortingEnabled(__sortingEnabled)
        self.Select_all.setText(_translate("Pondx", "Select all"))
        self.Clear_all.setText(_translate("Pondx", "Clear all"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Pondx = QtWidgets.QDialog()
    ui = Ui_Pondx()
    ui.setupUi(Pondx)
    Pondx.show()
    sys.exit(app.exec_())

