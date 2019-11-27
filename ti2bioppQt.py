# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ti2bioppQt.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(724, 659)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/newPrefix/RNA_2_Icon.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.listWidget = QtWidgets.QListWidget(self.centralwidget)
        self.listWidget.setObjectName("listWidget")
        self.verticalLayout_5.addWidget(self.listWidget)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.ShowFCM = QtWidgets.QPushButton(self.centralwidget)
        self.ShowFCM.setObjectName("ShowFCM")
        self.horizontalLayout.addWidget(self.ShowFCM)
        self.RemoveSeq = QtWidgets.QPushButton(self.centralwidget)
        self.RemoveSeq.setObjectName("RemoveSeq")
        self.horizontalLayout.addWidget(self.RemoveSeq)
        self.FastaSave = QtWidgets.QPushButton(self.centralwidget)
        self.FastaSave.setObjectName("FastaSave")
        self.horizontalLayout.addWidget(self.FastaSave)
        self.ClearBox = QtWidgets.QPushButton(self.centralwidget)
        self.ClearBox.setObjectName("ClearBox")
        self.horizontalLayout.addWidget(self.ClearBox)
        self.verticalLayout_5.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.Inputname = QtWidgets.QLineEdit(self.centralwidget)
        self.Inputname.setObjectName("Inputname")
        self.horizontalLayout_2.addWidget(self.Inputname)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        self.verticalLayout_5.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.inputSeq = QtWidgets.QLineEdit(self.centralwidget)
        self.inputSeq.setObjectName("inputSeq")
        self.horizontalLayout_3.addWidget(self.inputSeq)
        self.AddSeq = QtWidgets.QPushButton(self.centralwidget)
        self.AddSeq.setObjectName("AddSeq")
        self.horizontalLayout_3.addWidget(self.AddSeq)
        self.verticalLayout_5.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_4.addWidget(self.label_2)
        self.fileoutput = QtWidgets.QLineEdit(self.centralwidget)
        self.fileoutput.setObjectName("fileoutput")
        self.horizontalLayout_4.addWidget(self.fileoutput)
        self.DirChange = QtWidgets.QPushButton(self.centralwidget)
        self.DirChange.setObjectName("DirChange")
        self.horizontalLayout_4.addWidget(self.DirChange)
        self.verticalLayout_5.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.SelectoMode = QtWidgets.QGroupBox(self.centralwidget)
        self.SelectoMode.setObjectName("SelectoMode")
        self.layoutWidget = QtWidgets.QWidget(self.SelectoMode)
        self.layoutWidget.setGeometry(QtCore.QRect(34, 22, 156, 86))
        self.layoutWidget.setObjectName("layoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.layoutWidget)
        self.formLayout.setObjectName("formLayout")
        self.fcm = QtWidgets.QRadioButton(self.layoutWidget)
        self.fcm.setCheckable(True)
        self.fcm.setChecked(True)
        self.fcm.setObjectName("fcm")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.SpanningRole, self.fcm)
        self.nandy = QtWidgets.QRadioButton(self.layoutWidget)
        self.nandy.setObjectName("nandy")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.nandy)
        self.linealsort = QtWidgets.QRadioButton(self.layoutWidget)
        self.linealsort.setObjectName("linealsort")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.SpanningRole, self.linealsort)
        self.horizontalLayout_5.addWidget(self.SelectoMode)
        self.SelectorSort = QtWidgets.QGroupBox(self.centralwidget)
        self.SelectorSort.setObjectName("SelectorSort")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.SelectorSort)
        self.verticalLayout.setObjectName("verticalLayout")
        self.radio1 = QtWidgets.QRadioButton(self.SelectorSort)
        self.radio1.setChecked(True)
        self.radio1.setObjectName("radio1")
        self.verticalLayout.addWidget(self.radio1)
        self.radio2 = QtWidgets.QRadioButton(self.SelectorSort)
        self.radio2.setObjectName("radio2")
        self.verticalLayout.addWidget(self.radio2)
        self.radioCT = QtWidgets.QRadioButton(self.SelectorSort)
        self.radioCT.setObjectName("radioCT")
        self.verticalLayout.addWidget(self.radioCT)
        self.horizontalLayout_5.addWidget(self.SelectorSort)
        self.SetCodeway = QtWidgets.QGroupBox(self.centralwidget)
        self.SetCodeway.setObjectName("SetCodeway")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.SetCodeway)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.radio1_2 = QtWidgets.QRadioButton(self.SetCodeway)
        self.radio1_2.setChecked(True)
        self.radio1_2.setObjectName("radio1_2")
        self.verticalLayout_2.addWidget(self.radio1_2)
        self.radio2_2 = QtWidgets.QRadioButton(self.SetCodeway)
        self.radio2_2.setObjectName("radio2_2")
        self.verticalLayout_2.addWidget(self.radio2_2)
        self.horizontalLayout_5.addWidget(self.SetCodeway)
        self.verticalLayout_5.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.SetWeight = QtWidgets.QGroupBox(self.centralwidget)
        self.SetWeight.setObjectName("SetWeight")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.SetWeight)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.radioButton_5 = QtWidgets.QRadioButton(self.SetWeight)
        self.radioButton_5.setChecked(True)
        self.radioButton_5.setObjectName("radioButton_5")
        self.verticalLayout_3.addWidget(self.radioButton_5)
        self.radioButton_6 = QtWidgets.QRadioButton(self.SetWeight)
        self.radioButton_6.setObjectName("radioButton_6")
        self.verticalLayout_3.addWidget(self.radioButton_6)
        self.radioButton_7 = QtWidgets.QRadioButton(self.SetWeight)
        self.radioButton_7.setObjectName("radioButton_7")
        self.verticalLayout_3.addWidget(self.radioButton_7)
        self.radioButton_8 = QtWidgets.QRadioButton(self.SetWeight)
        self.radioButton_8.setObjectName("radioButton_8")
        self.verticalLayout_3.addWidget(self.radioButton_8)
        self.horizontalLayout_6.addWidget(self.SetWeight)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem1)
        self.morepond = QtWidgets.QPushButton(self.centralwidget)
        self.morepond.setObjectName("morepond")
        self.verticalLayout_4.addWidget(self.morepond)
        self.otherTIs = QtWidgets.QPushButton(self.centralwidget)
        self.otherTIs.setObjectName("otherTIs")
        self.verticalLayout_4.addWidget(self.otherTIs)
        self.gotoCalc = QtWidgets.QPushButton(self.centralwidget)
        self.gotoCalc.setObjectName("gotoCalc")
        self.verticalLayout_4.addWidget(self.gotoCalc)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem2)
        self.horizontalLayout_6.addLayout(self.verticalLayout_4)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.importRNA = QtWidgets.QPushButton(self.centralwidget)
        self.importRNA.setObjectName("importRNA")
        self.gridLayout.addWidget(self.importRNA, 4, 0, 1, 1)
        self.Cancel = QtWidgets.QPushButton(self.centralwidget)
        self.Cancel.setObjectName("Cancel")
        self.gridLayout.addWidget(self.Cancel, 4, 1, 1, 1)
        self.import_AAs = QtWidgets.QPushButton(self.centralwidget)
        self.import_AAs.setObjectName("import_AAs")
        self.gridLayout.addWidget(self.import_AAs, 3, 0, 1, 1)
        self.pdbimport = QtWidgets.QPushButton(self.centralwidget)
        self.pdbimport.setObjectName("pdbimport")
        self.gridLayout.addWidget(self.pdbimport, 3, 1, 1, 1)
        self.importSeqNA = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.importSeqNA.sizePolicy().hasHeightForWidth())
        self.importSeqNA.setSizePolicy(sizePolicy)
        self.importSeqNA.setObjectName("importSeqNA")
        self.gridLayout.addWidget(self.importSeqNA, 1, 0, 1, 1)
        self.importCT = QtWidgets.QPushButton(self.centralwidget)
        self.importCT.setObjectName("importCT")
        self.gridLayout.addWidget(self.importCT, 1, 1, 1, 1)
        self.horizontalLayout_6.addLayout(self.gridLayout)
        self.verticalLayout_5.addLayout(self.horizontalLayout_6)
        self.progressBar = QtWidgets.QProgressBar(self.centralwidget)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setTextVisible(False)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout_5.addWidget(self.progressBar)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 724, 25))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(MainWindow)
        self.toolBar.setFloatable(True)
        self.toolBar.setObjectName("toolBar")
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionPDB = QtWidgets.QAction(MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/newPrefix/mydrawing.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPDB.setIcon(icon1)
        self.actionPDB.setObjectName("actionPDB")
        self.actionDNA_import = QtWidgets.QAction(MainWindow)
        self.actionDNA_import.setCheckable(True)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/newPrefix/DNA1.PNG"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionDNA_import.setIcon(icon2)
        self.actionDNA_import.setWhatsThis("")
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.actionDNA_import.setFont(font)
        self.actionDNA_import.setIconVisibleInMenu(True)
        self.actionDNA_import.setObjectName("actionDNA_import")
        self.actionProtein_import = QtWidgets.QAction(MainWindow)
        self.actionProtein_import.setCheckable(True)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/newPrefix/protein.PNG"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionProtein_import.setIcon(icon3)
        self.actionProtein_import.setObjectName("actionProtein_import")
        self.actionClose = QtWidgets.QAction(MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/newPrefix/Exit_icon.PNG"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionClose.setIcon(icon4)
        self.actionClose.setObjectName("actionClose")
        self.actionCalculate = QtWidgets.QAction(MainWindow)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(":/newPrefix/compute_icon.PNG"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionCalculate.setIcon(icon5)
        self.actionCalculate.setObjectName("actionCalculate")
        self.toolBar.addAction(self.actionDNA_import)
        self.toolBar.addAction(self.actionProtein_import)
        self.toolBar.addAction(self.actionPDB)
        self.toolBar.addAction(self.actionCalculate)
        self.toolBar.addAction(self.actionClose)

        self.retranslateUi(MainWindow)
        self.Cancel.clicked.connect(MainWindow.close)
        self.actionClose.triggered.connect(MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "TI2BioP 2.1"))
        self.ShowFCM.setText(_translate("MainWindow", "Show"))
        self.RemoveSeq.setText(_translate("MainWindow", "Remove selected"))
        self.FastaSave.setText(_translate("MainWindow", "Save"))
        self.ClearBox.setText(_translate("MainWindow", "Clear box"))
        self.Inputname.setText(_translate("MainWindow", "noname_0"))
        self.label.setText(_translate("MainWindow", "Sequence name"))
        self.AddSeq.setText(_translate("MainWindow", "Add"))
        self.label_2.setText(_translate("MainWindow", "Output filename"))
        self.fileoutput.setText(_translate("MainWindow", "noname.txt"))
        self.DirChange.setText(_translate("MainWindow", "..Dir.."))
        self.SelectoMode.setTitle(_translate("MainWindow", "Apply"))
        self.fcm.setText(_translate("MainWindow", "Four color maps"))
        self.nandy.setText(_translate("MainWindow", "Nandy Graph"))
        self.linealsort.setText(_translate("MainWindow", "lineal"))
        self.SelectorSort.setTitle(_translate("MainWindow", "Proccess as:"))
        self.radio1.setText(_translate("MainWindow", "DNA/RNA"))
        self.radio2.setText(_translate("MainWindow", "Amino acids"))
        self.radioCT.setText(_translate("MainWindow", "CT with Secondary Str."))
        self.SetCodeway.setTitle(_translate("MainWindow", "Amino acids code"))
        self.radio1_2.setText(_translate("MainWindow", "One Letter (A, C, F ...)"))
        self.radio2_2.setText(_translate("MainWindow", "Three Letters (Ala, Cys, Phe...)"))
        self.SetWeight.setTitle(_translate("MainWindow", "Select Weight"))
        self.radioButton_5.setText(_translate("MainWindow", "sum-Amber95"))
        self.radioButton_6.setText(_translate("MainWindow", "sum-FCPAC free"))
        self.radioButton_7.setText(_translate("MainWindow", "sum-FCPAC buried"))
        self.radioButton_8.setText(_translate("MainWindow", "ISA"))
        self.morepond.setText(_translate("MainWindow", "Ponderations...."))
        self.otherTIs.setText(_translate("MainWindow", "Other TIs"))
        self.gotoCalc.setText(_translate("MainWindow", "TIs Calculation"))
        self.importRNA.setText(_translate("MainWindow", "Import RNA xfasta"))
        self.Cancel.setText(_translate("MainWindow", "&Close"))
        self.Cancel.setShortcut(_translate("MainWindow", "Alt+C"))
        self.import_AAs.setText(_translate("MainWindow", "Import AA Seq"))
        self.pdbimport.setText(_translate("MainWindow", "Pdbs"))
        self.importSeqNA.setText(_translate("MainWindow", "Import DNA Seq"))
        self.importCT.setText(_translate("MainWindow", "Import CT"))
        self.toolBar.setWindowTitle(_translate("MainWindow", "TI2BioP 2.0"))
        self.actionPDB.setText(_translate("MainWindow", "PDB"))
        self.actionPDB.setStatusTip(_translate("MainWindow", "Import pdb files"))
        self.actionPDB.setShortcut(_translate("MainWindow", "Ctrl+V"))
        self.actionDNA_import.setText(_translate("MainWindow", "DNA import"))
        self.actionDNA_import.setIconText(_translate("MainWindow", "DNA import", "DNA import"))
        self.actionDNA_import.setToolTip(_translate("MainWindow", "DNA Import"))
        self.actionDNA_import.setStatusTip(_translate("MainWindow", "Import DNA/RNA sequences"))
        self.actionDNA_import.setShortcut(_translate("MainWindow", "Alt+D"))
        self.actionProtein_import.setText(_translate("MainWindow", "Protein import"))
        self.actionProtein_import.setToolTip(_translate("MainWindow", "Protein import"))
        self.actionProtein_import.setStatusTip(_translate("MainWindow", "Import Protein Sequences"))
        self.actionProtein_import.setShortcut(_translate("MainWindow", "Alt+P"))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionClose.setToolTip(_translate("MainWindow", "Close Application"))
        self.actionClose.setStatusTip(_translate("MainWindow", "Close application"))
        self.actionCalculate.setText(_translate("MainWindow", "Calculate"))
        self.actionCalculate.setToolTip(_translate("MainWindow", "Calculations "))
        self.actionCalculate.setStatusTip(_translate("MainWindow", "Perform Calculations"))
        self.actionCalculate.setShortcut(_translate("MainWindow", "Alt+Shift+C"))

import reso_rc

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
