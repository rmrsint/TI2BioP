# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ponderation_NA.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_PondRNA(object):
    def setupUi(self, PondRNA):
        PondRNA.setObjectName("PondRNA")
        PondRNA.resize(568, 279)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/newPrefix/RNA_2_Icon.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        PondRNA.setWindowIcon(icon)
        self.verticalLayout = QtWidgets.QVBoxLayout(PondRNA)
        self.verticalLayout.setObjectName("verticalLayout")
        self.listWidget = QtWidgets.QListWidget(PondRNA)
        self.listWidget.setObjectName("listWidget")
        item = QtWidgets.QListWidgetItem()
        item.setCheckState(QtCore.Qt.Checked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
        item.setCheckState(QtCore.Qt.Unchecked)
        self.listWidget.addItem(item)
        item = QtWidgets.QListWidgetItem()
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
        self.Select_all = QtWidgets.QPushButton(PondRNA)
        self.Select_all.setObjectName("Select_all")
        self.horizontalLayout.addWidget(self.Select_all)
        self.Clear_all = QtWidgets.QPushButton(PondRNA)
        self.Clear_all.setObjectName("Clear_all")
        self.horizontalLayout.addWidget(self.Clear_all)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        self.buttonBox = QtWidgets.QDialogButtonBox(PondRNA)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.horizontalLayout_2.addWidget(self.buttonBox)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(PondRNA)
        self.buttonBox.accepted.connect(PondRNA.accept)
        self.buttonBox.rejected.connect(PondRNA.reject)
        QtCore.QMetaObject.connectSlotsByName(PondRNA)

    def retranslateUi(self, PondRNA):
        _translate = QtCore.QCoreApplication.translate
        PondRNA.setWindowTitle(_translate("PondRNA", "Ponderations"))
        __sortingEnabled = self.listWidget.isSortingEnabled()
        self.listWidget.setSortingEnabled(False)
        item = self.listWidget.item(0)
        item.setText(_translate("PondRNA", "sum-Amber95"))
        item = self.listWidget.item(1)
        item.setText(_translate("PondRNA", "first oscillator strength value ( singlet excitation energies)"))
        item = self.listWidget.item(2)
        item.setText(_translate("PondRNA", "second oscillator strength value"))
        item = self.listWidget.item(3)
        item.setText(_translate("PondRNA", "molar absorption coefficient at 260 nm and PH = 7.0,"))
        item = self.listWidget.item(4)
        item.setText(_translate("PondRNA", "first (ΔE1) single excitation energies in eV"))
        item = self.listWidget.item(5)
        item.setText(_translate("PondRNA", "second (ΔE2) single excitation energies in eV"))
        self.listWidget.setSortingEnabled(__sortingEnabled)
        self.Select_all.setText(_translate("PondRNA", "Select all"))
        self.Clear_all.setText(_translate("PondRNA", "Clear all"))

import reso_rc

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    PondRNA = QtWidgets.QDialog()
    ui = Ui_PondRNA()
    ui.setupUi(PondRNA)
    PondRNA.show()
    sys.exit(app.exec_())

