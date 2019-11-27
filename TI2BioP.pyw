#!/usr/bin/python3
__author__ = 'reymolina'

global app
import TI2BioP
import Bio
import TI2BioP_nandy
import openct as ct
import openxfasta as xfasta
import sqlalchemy as sq
import csv
import time
import re, os
import numpy as np
from ti2bioppQt import *
from oTIsDialog import *
from Bio.SeqUtils import seq1

from PyQt5 import QtGui, QtCore, QtWidgets

from ponddialog import *
from pondNAdialog import *


def checkEqual(lst):
    return lst[1:] == lst[:-1]


def fixSequences(sequences):
    NewSequences = []
    for x in sequences:
        NewSequences.append([x['name'], x['seq']])
    return NewSequences


def getprot(filename):
    pattern = r'^ATOM *\d+  CA  (?P<atom>\w{3})'
    regexp = re.compile(pattern)

    pdb = open(filename, 'r')

    aminoacid = []

    for linea in pdb:
        result = regexp.search(linea)
        if result:
            # print (result.group('atom'),'--->',seq1(result.group('atom')) )
            aminoacid.append(seq1(result.group('atom')))

    return ''.join(aminoacid)


def checkprot(s, codesort):
    '''
    :param s:
    :param codesort :
     0 one letter code for amino acids
     1 three letters code
    :return boolean:
    '''

    marker = True
    if codesort == 0:
        for i in range(len(s)):
            if s[i] in ['O', 'U', 'J']:
                marker = False
                break
    elif codesort == 1:
        s1 = seq1(s)
        if s1.count('X') != 0:
            marker = False
        else:
            marker = True

    return marker


def checkNucAcids(s):
    marker = True
    for i in range(len(s)):
        if s[i] not in ['A', 'T', 'U', 'G', 'C']:
            marker = False
            break
    return marker


class TaskLoadfile(QtCore.QThread):
    taskFinished = QtCore.pyqtSignal()
    total = QtCore.pyqtSignal(object)

    def __init__(self, parent, *args, **kw):
        super(TaskLoadfile, self).__init__(parent)
        self.myInit(*args, **kw)
        self.results = []

    def myInit(self, filenames):
        self.filenames = filenames

    def run(self):
        time.sleep(3)
        self.taskFinished.emit()


class TaskHeavy(QtCore.QThread):
    taskFinished = QtCore.pyqtSignal()
    total = QtCore.pyqtSignal(object)

    def __init__(self, parent, *args, **kw):
        super(TaskHeavy, self).__init__(parent)
        self.myInit(*args, **kw)
        self.results = []

    def myInit(self, sequences, ponderations, consider_as):
        self.sequences = sequences
        self.ponderations = ponderations
        self.consider_as = consider_as

    def setorientation(self, orientation):
        self.orientation = orientation

    def run(self):
        time.sleep(3)
        self.taskFinished.emit()


class calculation_1(QtCore.QThread):
    total = QtCore.pyqtSignal(object)
    update = QtCore.pyqtSignal()

    def __init__(self, parent, *args, **kw):
        super(calculation_1, self).__init__(parent)
        self.myInit(*args, **kw)
        self.results = []

    def myInit(self, sequences, ponderations, consider_as):
        self.sequences = sequences
        self.ponderations = ponderations
        self.consider_as = consider_as

    def run(self):
        self.total.emit(self.sequences.__len__())
        for i in range(len(self.sequences)):
            macromolecule = TI2BioP.MacroMol(self.sequences[i], self.ponderations, self.consider_as[i], True)
            self.results.append(macromolecule.result)
            self.update.emit()


class progress1(QtWidgets.QProgressBar):

    def __init__(self, parent=None):
        super(progress1, self).__init__(parent)
        # Set up the user interface from Designer.
        self.setValue(0)

    def prepareready(self, sequences, ponderations, consider_as):
        self.sequences = sequences
        self.ponderations = ponderations
        self.consider_as = consider_as

        self.thread = calculation_1(self, self.sequences, self.ponderations, self.consider_as)
        self.thread.total.connect(self.setMaximum)
        self.thread.update.connect(self.update)
        self.thread.finished.connect(self.close)

        self.n = 0
        self.thread.start()

    def update(self):
        self.n += 1
        print (self.n)
        self.setValue(self.n)


class MyInterface(QtWidgets.QMainWindow):
    signal = QtCore.pyqtSignal(object)


    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ponderations = [0]
        self.filenamex = 'noname'
        self.ui.AddSeq.clicked.connect(self.addItem)
        self.ui.ClearBox.clicked.connect(self.clearBox)
        self.ui.RemoveSeq.clicked.connect(self.delItem)
        self.ui.morepond.clicked.connect(self.morepond)
        self.ui.gotoCalc.clicked.connect(self.preprocessing)
        self.ui.actionCalculate.triggered.connect(self.preprocessing)
        self.ui.inputSeq.returnPressed.connect(self.addItem)
        self.ui.importSeqNA.clicked.connect(self.importNA)
        self.ui.importSeqNA.clicked.connect(self.importNA)
        self.ui.importRNA.clicked.connect(self.importRNAx)
        self.ui.actionDNA_import.triggered.connect(self.importNA)
        self.ui.import_AAs.clicked.connect(self.importAAs)
        self.ui.actionProtein_import.triggered.connect(self.importAAs)
        self.ui.DirChange.clicked.connect(self.changedir)
        self.ui.ShowFCM.clicked.connect(self.showspiral)
        self.ui.actionPDB.triggered.connect(self.importPdbs)
        self.ui.radioButton_5.clicked.connect(self.changepond)
        self.ui.radioButton_6.clicked.connect(self.changepond)
        self.ui.radioButton_7.clicked.connect(self.changepond)
        self.ui.radioButton_8.clicked.connect(self.changepond)
        self.ui.pdbimport.clicked.connect(self.importPdbs)
        self.ui.importCT.clicked.connect(self.importCT)
        self.ui.FastaSave.clicked.connect(self.savefasta)
        self.ui.otherTIs.clicked.connect(self.otherTIs)



    def update(self):
        self.n += 1
        self.ui.progressBar.setValue(self.n)

    def reciveresult(self, total):
        print(total)

    def changedir(self):
        dirx = QtWidgets.QFileDialog.getExistingDirectory(None, 'Select a folder:', '.', QtWidgets.QFileDialog.ShowDirsOnly)
        temp = str(self.ui.fileoutput.text())
        if '/' in temp:
            temp = temp.split('/')[-1]

        self.ui.fileoutput.setText(dirx + '/' + temp)

    def stoploading(self):
        # stop the pulsation
        for x in self.results:
            self.ui.listWidget.addItem(x)
        self.results = []

        self.ui.progressBar.setRange(0, 1)
        self.ui.statusbar.showMessage("the files were loaded and proccessed.")

    def stopcalculation(self):
        # stop the pulsation
        self.ui.progressBar.setRange(0, 1)
        self.ui.statusbar.showMessage("Calculations were performed. Please check your output file")
        QtWidgets.QMessageBox.information(self, "Calculation Successful", "The calculations were performed ...")

    def changepond(self):
        sender = self.sender()
        self.statusBar().showMessage(sender.text() + ' ponderation was selected')
        # print sender.objectName()
        podn = int(str(sender.objectName()).split('_')[1]) - 5
        self.ponderations = [podn]
        # print self.ponderations

    def calcGlobal(self):
        for i in range(len(self.sequences)):
            macromolecule = TI2BioP.MacroMol(self.sequences[i], self.ponderations, self.consider_as[i], True)
            self.results.append(macromolecule.result)
        self.signal.emit(self.results)
        self.saveoutput()
        self.myLongTask.taskFinished.emit()

    def calcGlobal_all(self):
        if len(self.sequences[0]) > 2:
            self.sequences = fixSequences(self.sequences)

        self.results = TI2BioP.process_all(self.sequences, self.ponderations, self.consider_as[0])
        self.signal.emit(self.results)
        self.saveoutput()
        self.myLongTask.taskFinished.emit()

    def loadprocessfile(self):
        self.results = []
        # input and proccess files
        fileinput = self.filenames
        for i in fileinput:
            basename = os.path.basename(str(i))
            seqdata = xfasta.allfromxfasta(str(i))
            for j in range(len(seqdata.xfasta)):
                self.results.append('>' + basename + '|' + seqdata.xfasta[j][0])
                self.results.append(seqdata.xfasta[j][1] + '|Nucleic acids|' + str(seqdata.xfasta[j][2]))

        self.signal.emit(self.results)
        self.myLongTask.taskFinished.emit()

    def calcGlobalNandy_all_lineal(self):
        if len(self.sequences[0]) > 2:
            self.sequences = fixSequences(self.sequences)

        if self.consider_as[0] == 2:
            self.results = TI2BioP_nandy.process_all(self.sequences, TI2BioP_nandy.orienta_seq_lineal,
                                                     self.ponderations, self.consider_as[0])
        else:
            self.results = TI2BioP_nandy.process_all(self.sequences, TI2BioP_nandy.orienta_amino_lineal,
                                                     self.ponderations, self.consider_as[0])
        self.signal.emit(self.results)
        self.saveoutput()
        self.myLongTask.taskFinished.emit()

    def calcGlobalNandy_all(self):
        if len(self.sequences[0]) > 2:
            self.sequences = fixSequences(self.sequences)

        if self.consider_as[0] == 2:
            self.results = TI2BioP_nandy.process_all(self.sequences, TI2BioP_nandy.orientaseq, self.ponderations,
                                                     self.consider_as[0])
        else:
            self.results = TI2BioP_nandy.process_all(self.sequences, TI2BioP_nandy.orientAmino, self.ponderations,
                                                     self.consider_as[0])

        self.signal.emit(self.results)
        self.saveoutput()
        self.savetosql()

        self.myLongTask.taskFinished.emit()

    def calcGlobalNandy_one_by_one(self):
        if len(self.sequences[0]) > 2:
            self.sequences = fixSequences(self.sequences)

        for i in range(len(self.sequences)):
            macromolecule = TI2BioP_nandy.MacroMol(self.sequences[i], self.ponderations, self.consider_as[i], True)
            self.results.append(macromolecule.result)

        self.signal.emit(self.results)
        self.saveoutput()
        self.myLongTask.taskFinished.emit()

    def calcGlobalNandy_1_1_lineal(self):
        if len(self.sequences[0]) > 2:
            self.sequences = fixSequences(self.sequences)

        for i in range(len(self.sequences)):
            if self.consider_as[i] == 2:
                macromolecule = TI2BioP_nandy.MacroMol(self.sequences[i], TI2BioP_nandy.orienta_seq_lineal,
                                                       self.ponderations, 1, False, True)
            else:
                macromolecule = TI2BioP_nandy.MacroMol(self.sequences[i], TI2BioP_nandy.orienta_amino_lineal,
                                                       self.ponderations, 0, False, True)
            self.results.append(macromolecule.result)

        self.signal.emit(self.results)
        self.saveoutput()
        self.myLongTask.taskFinished.emit()

    def calc_all_CT(self):
        self.results = ct.process_all(self.sequences, self.ponderations)
        self.signal.emit(self.results)

        self.saveoutput()
        self.myLongTask.taskFinished.emit()

    def center(self):
        frameGm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerPoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())

    def addItem(self):
        process_as = ''
        s = str(self.ui.inputSeq.text()).upper()
        c = s.split()
        stocheck = ''.join(c)
        additem_x = False

        if self.ui.radio1.isChecked():
            additem_x = checkNucAcids(stocheck)
            process_as = 'Nucleic acids'

        if self.ui.radio2.isChecked():
            process_as = 'Protein'
            if self.ui.radio1_2.isChecked():
                additem_x = checkprot(stocheck, 0)
            elif self.ui.radio2_2.isChecked():
                additem_x = checkprot(stocheck, 1)
                if additem_x:
                    stocheck = seq1(stocheck)

        if additem_x:
            if len(self.ui.Inputname.text()) == 0:
                self.ui.Inputname.setText('seq_' + str(self.ui.listWidget.count() // 2))

            if '>' not in self.ui.Inputname.text():
                self.ui.listWidget.addItem('>' + self.ui.Inputname.text())
            else:
                self.ui.listWidget.addItem(self.ui.Inputname.text())

            qtfield = stocheck + '|' + process_as

            self.ui.listWidget.addItem(qtfield)
            # self.ui.listWidget.addItem(QtCore.QString('Process as: ' + process_as))
            self.ui.Inputname.setText('')
            self.ui.inputSeq.setText('')
            self.ui.Inputname.setFocus()
            self.ui.inputSeq.setFocus()
        else:
            reply = QtWidgets.QMessageBox.question(self, 'Message  Sequence Error !!!',
                                               'Do you want to edit the input sequence ?',
                                               QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if reply == QtWidgets.QMessageBox.No:
                self.ui.inputSeq.setText('')
                self.ui.inputSeq.setFocus()

    def delItem(self):
        row = self.ui.listWidget.currentRow()
        if row % 2 == 0:
            self.ui.listWidget.takeItem(row)
            self.ui.listWidget.takeItem(row + 1)
        else:
            self.ui.listWidget.takeItem(row)
            self.ui.listWidget.takeItem(row - 1)

    def clearBox(self):
        self.ui.listWidget.clear()

    def morepond(self):
        if self.ui.radio1.isChecked() or self.ui.radioCT.isChecked():
            PonderationDialog = DialogPondRNA()
            if PonderationDialog.exec_():
                self.ponderations = PonderationDialog.getvalues()
            print(self.ponderations)
        else:
            PonderationDialog = DialogPond()
            PonderationDialog.Addponderations()
            if PonderationDialog.exec_():
                self.ponderations = PonderationDialog.getvalues()

    def showspiral(self):
        row = self.ui.listWidget.currentRow()
        if row != -1:
            if row % 2 == 0:
                namex = str(self.ui.listWidget.item(row).text())
                item = self.ui.listWidget.item(row + 1)
            else:
                item = self.ui.listWidget.item(row)
                namex = str(self.ui.listWidget.item(row - 1).text())

            s = str(item.text())
            c = s.split('|')

            # print c[1]

            if c[1] == 'Nucleic acids':
                if 'U' in c[0]:
                    see_as = 3
                else:
                    see_as = 2
            else:
                see_as = 0

            # print ('----> ',see_as)
            if self.ui.radio1.isChecked():
                self.ponderations = [0]

            if self.ui.fcm.isChecked():
                macromolecule = TI2BioP.MacroMol([namex, c[0]], [], see_as, False)
                # To show four colors map representation
                macromolecule.showspiral2()

            if self.ui.nandy.isChecked():
                if self.ui.radio1.isChecked():
                    macromolecule = TI2BioP_nandy.MacroMol([namex, c[0]], TI2BioP_nandy.orientaseq, [], 2, True, False)
                else:
                    if self.ui.radio1_2.isChecked():
                        macromolecule = TI2BioP_nandy.MacroMol([namex, c[0]], TI2BioP_nandy.orientaseq, [], 0, True,
                                                               False)
                    else:
                        macromolecule = TI2BioP_nandy.MacroMol([namex, c[0]], TI2BioP_nandy.orientaseq, [], 1, True,
                                                               False)

            if self.ui.linealsort.isChecked():
                if self.ui.radio1.isChecked():
                    macromolecule = TI2BioP_nandy.MacroMol([namex, c[0]], TI2BioP_nandy.orienta_seq_lineal, [], 0, True,
                                                           False)

                if self.ui.radio2.isChecked():
                    if self.ui.radio1_2:
                        macromolecule = TI2BioP_nandy.MacroMol([namex, c[0]], TI2BioP_nandy.orienta_amino_lineal, [], 0,
                                                               True, False)
                    if self.ui.radio2_2.isChecked():
                        macromolecule = TI2BioP_nandy.MacroMol([namex, c[0]], TI2BioP_nandy.orienta_seq_lineal, [], 1,
                                                               True, False)

    def preprocessing(self):
        calculation = False
        self.names = []
        self.consider_as = []
        self.sequences = []

        if self.ui.listWidget.count() == 0:
            self.ui.statusbar.showMessage("No calculations could be performed | check your inputs")
            QtWidgets.QMessageBox.critical(self, "Invalid request", "You must input at the least one sequence to calculate",
                                       QtWidgets.QMessageBox.Ok)

        else:
            calculation = True
            self.ui.statusbar.showMessage("Processing calculations ...")

            for i in range(self.ui.listWidget.count()):
                item = self.ui.listWidget.item(i)
                s = str(item.text())
                if '>' in s:
                    names = str(s.split('>')[1])
                else:
                    c = s.split('|')
                    if len(c) == 2:
                        self.sequences.append([names, str(c[0])])
                    else:
                        self.sequences.append({'name': names, 'seq': c[0], 'connection': c[2]})

                    if c[1] == 'Nucleic acids':
                        if 'U' in c[0]:
                            self.consider_as.append(2)
                        else:
                            self.consider_as.append(2)
                    else:
                        self.consider_as.append(0)

        self.ui.progressBar.setRange(0, 0)
        self.ui.progressBar.setValue(0)
        if calculation:
            self.results = []

            if self.ui.radioCT.isChecked():
                print(self.ponderations, '<----->')
                self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                self.myLongTask.taskFinished.connect(self.stopcalculation)
                self.myLongTask.run = self.calc_all_CT
                self.myLongTask.start()
            else:
                if self.ui.fcm.isChecked():
                    if checkEqual(list(self.consider_as)):
                        self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                        self.myLongTask.taskFinished.connect(self.stopcalculation)
                        self.myLongTask.run = self.calcGlobal_all
                        self.myLongTask.start()
                        # self.results = TI2BioP.process_all(self.sequences,self.ponderations, self.consider_as[0])
                        # self.ui.statusbar.showMessage("Calculations were performed. Please check your output file")
                    else:
                        self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                        self.myLongTask.taskFinished.connect(self.stopcalculation)
                        self.myLongTask.run = self.calcGlobal
                        self.myLongTask.start()

                if self.ui.nandy.isChecked():

                    if checkEqual(list(self.consider_as)):
                        self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                        self.myLongTask.taskFinished.connect(self.stopcalculation)
                        self.myLongTask.run = self.calcGlobalNandy_all
                        self.myLongTask.start()
                    else:
                        self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                        self.myLongTask.taskFinished.connect(self.stopcalculation)
                        self.myLongTask.run = self.calcGlobalNandy_one_by_one
                        self.myLongTask.start()
                        '''
                        for i in range(len(self.sequences)):
                            app.processEvents()
                            #hay que arreglar
                            macromolecule = TI2BioP_nandy.MacroMol(self.sequences[i],self.ponderations, self.consider_as[i], True)
                            self.results.append(macromolecule.result)
                            self.ui.progressBar.setValue(i * 100.0 / len(self.sequences))
                        self.ui.statusbar.showMessage("Calculations were performed | check your output file")
                        '''

                if self.ui.linealsort.isChecked():

                    if checkEqual(list(self.consider_as)):
                        self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                        self.myLongTask.taskFinished.connect(self.stopcalculation)
                        self.myLongTask.run = self.calcGlobalNandy_all_lineal
                        self.myLongTask.start()
                        '''
                        if self.consider_as[0] == 2:
                            self.results = TI2BioP_nandy.process_all(self.sequences,TI2BioP_nandy.orienta_seq_lineal,self.ponderations, self.consider_as[0])
                        else:
                            self.results = TI2BioP_nandy.process_all(self.sequences,TI2BioP_nandy.orienta_amino_lineal,self.ponderations, self.consider_as[0])
                        self.ui.statusbar.showMessage("Calculations were performed | check your output file")
                        '''
                    else:

                        self.myLongTask = TaskHeavy(None, self.sequences, self.ponderations, self.consider_as)
                        self.myLongTask.taskFinished.connect(self.stopcalculation)
                        self.myLongTask.run = self.calcGlobalNandy_1_1_lineal
                        self.myLongTask.start()

    def importRNAx(self):
        fileinput = QtWidgets.QFileDialog.getOpenFileNames(self, 'Choice data file...', '', 'XFasta file(*.xfasta)')
        self.filenames = fileinput[0]
        self.ui.progressBar.setRange(0, 0)
        self.ui.progressBar.setValue(0)

        if fileinput[0]:
            self.myLongTask = TaskLoadfile(None, self.filenames)
            self.myLongTask.taskFinished.connect(self.stoploading)
            self.myLongTask.run = self.loadprocessfile
            self.myLongTask.start()

        self.ui.radio1.setChecked(True)
        self.ui.radioCT.setChecked(True)

        #    self.ui.statusbar.showMessage("All xfasta were processed ...")
        # else:
        #    self.ui.statusbar.showMessage("No xfasta file  to process ...")

    def importNA(self):
        fileinput = QtWidgets.QFileDialog.getOpenFileName(self, 'Choice data file...', '', 'Fasta file(*.fasta)')

        if fileinput[0]:
            f = open(str(fileinput[0]), 'r')
            lines = f.readlines()
            cad = ""
            for x in lines:
                if x != '\n':
                    if ('>' in x) or ('[' in x) or (';' in x):
                        if len(cad) != 0:
                            self.ui.listWidget.addItem(seqname)
                            self.ui.listWidget.addItem(cad + '|Nucleic acids')
                            seqname = x.strip()
                            cad = ""
                        else:
                            seqname = x.strip()
                    else:
                        cad = cad + x.strip()
            self.ui.listWidget.addItem(seqname)
            self.ui.listWidget.addItem(cad + '|Nucleic acids')
        self.ui.radio1.setChecked(True)

    def importPdbs(self):

        fileinput = QtWidgets.QFileDialog.getOpenFileNames(self, 'Choice data file...', '', 'pdb files(*.pdb *.ent)')

        if fileinput[0]:
            for i in fileinput[0]:
                basename = os.path.basename(str(i))
                self.ui.statusbar.showMessage("Processing :  " + basename)
                proteinseq = getprot(i)
                self.ui.listWidget.addItem('>' + basename)
                self.ui.listWidget.addItem(proteinseq + '|Protein')

            self.ui.radio2.setChecked(True)

            self.ui.statusbar.showMessage("All pdbs were processed ...")
        else:
            self.ui.statusbar.showMessage("No pdbs to process ...")

    def importCT(self):

        fileinput = QtWidgets.QFileDialog.getOpenFileNames(self, 'Choice data file...', '', 'ct files(*.ct)')

        if fileinput[0]:
            for i in fileinput[0]:
                print(i)
                basename = os.path.basename(str(i))
                self.ui.statusbar.showMessage("Processing :  " + basename)
                seq = ct.allfromctrecord(i)
                for j in range(len(seq.ct)):
                    self.ui.listWidget.addItem('>' + basename + '|' + seq.ct[j][0])
                    self.ui.listWidget.addItem(seq.ct[j][1] + '|Nucleic acids|' + str(seq.ct[j][2]))

            self.ui.radio1.setChecked(True)
            self.ui.radioCT.setChecked(True)

            self.ui.statusbar.showMessage("All ct were processed ...")
        else:
            self.ui.statusbar.showMessage("No ct file  to process ...")

    def importAAs(self):
        patterns = ['>', '[', ';']

        fileinput = QtWidgets.QFileDialog.getOpenFileName(self, 'Choice data file...', '',
                                                      'Fasta file(*.fasta);;Text Files (*.txt);;All Files (*)',
                                                      'small data base (*.txt)')
        if fileinput[0]:
            self.filenamex = fileinput[0]
            f = open(str(fileinput[0]), 'r')
            lines = f.readlines()
            cad = ""
            for x in lines:
                if x != '\n':
                    if ('>' in x) or ('[' in x) or (';' in x):
                        if len(cad) != 0:
                            self.ui.listWidget.addItem(seqname)
                            self.ui.listWidget.addItem(cad + '|Protein')
                            seqname = x.strip()
                            cad = ""
                        else:
                            seqname = x.strip()
                    else:
                        cad = cad + x.strip()
            self.ui.listWidget.addItem(seqname)
            self.ui.listWidget.addItem(cad + '|Protein')
        self.ui.radio2.setChecked(True)

    def savefasta(self):
        if self.ui.listWidget.count() != 0:
            fileforsave = QtWidgets.QFileDialog.getSaveFileName(self, 'Save as ...', QtCore.QDir.currentPath(),
                                                           'Fasta files (*.fasta);;All files(*)')
            filetosave = str(fileforsave[0])

            if filetosave:
                if not '.' in filetosave:
                    filetosave = filetosave + '.fasta'

                f = open(filetosave, 'w')
                for i in range(self.ui.listWidget.count()):
                    item = self.ui.listWidget.item(i)
                    s = str(item.text())
                    if '>' in s:
                        f.write(s + '\n')

                    if '|' in s:
                        seqs = str(s.split('|')[0])
                        f.write(seqs + '\n')
                f.close()
            self.ui.statusbar.showMessage("All pbd files were saved...")
        else:
            self.ui.statusbar.showMessage("There are not sequences to save ...")

    def savetosql(self):

        # print (sq.__version__)

        if self.ui.radioCT.isChecked():
            suffix = '_sCT_'

            strk = "Name" + "\t" + "Edge Number" + "\t" + "Edge connnectivity Index" + "\t" + "Ponderation" + "\t"
            for i in range(1, 16):
                strk = strk + "Momment" + str(i) + "\t"

            strk = strk + '\n'
            name_out = str(self.ui.fileoutput.text()).split('.')
            namefinal = name_out[0] + suffix + '.db'

        else:

            if self.ui.fcm.isChecked():
                suffix = '_fcm_'
                strk = 'Name\tCluster Number\tTotal adjacency index\tPonderation\t'
            else:
                suffix = '_nandy_'
                strk = 'Name\tTotal adjacency index\tPonderation\t'

            if self.ui.linealsort.isChecked():
                suffix = '_lineal_'

            for i in range(1, 16):
                strk = strk + 'Moment ' + str(i) + '\t'

            # print(suffix)
            # print self.results
            # print len(self.results)
            # print self.ponderations
            # print len(self.results[0])

            if self.ponderations.__len__() != 1:
                seedName = str(self.ui.fileoutput.text()).split('.')
                for i in range(self.ponderations.__len__()):
                    name_out = seedName[0] + suffix + str(i) + '.db'
                    asalvar = []
                    for j in range(self.results.__len__()):

                        if len(self.results[j]) - self.ponderations.__len__() == 1:

                            asalvar.append([self.results[j][0]] + list(self.results[j][i + 1]))
                        elif len(self.results[j]) - self.ponderations.__len__() == 2:

                            asalvar.append([self.results[j][0]] + [self.results[j][1]] + list(self.results[j][i + 2]))

                    np.savetxt(name_out, asalvar, delimiter='\t', header=strk, fmt='%s')
            else:
                name_out = str(self.ui.fileoutput.text()).split('.')

                asalvar = []

                for i in range(self.results.__len__()):
                    if len(self.results[i]) == 2:
                        asalvar.append([self.results[i][0]] + list(self.results[i][1]))
                    elif len(self.results[i]) == 3:
                        asalvar.append([self.results[i][0]] + [self.results[i][1]] + list(self.results[i][2]))

                np.savetxt(name_out[0] + suffix + '.db', asalvar, delimiter='\t', header=strk, fmt='%s')

    def otherTIs(self):
        tpcalculus = DialogTIs()

        self.ui.statusbar.showMessage("processing calculations ...")
        xp = []
        consider_asX = []

        for i in range(self.ui.listWidget.count()):
            item = self.ui.listWidget.item(i)
            s = str(item.text())
            if '>' in s:
                names = str(s.split('>')[1])
            else:
                c = s.split('|')
                if len(c) == 2:
                    xp.append([names, str(c[0])])
                else:
                    xp.append({'name': names, 'seq': c[0], 'connection': c[2]})

                if c[1] == 'Nucleic acids':
                    if 'U' in c[0]:
                        consider_asX.append(2)
                    else:
                        consider_asX.append(2)
                else:
                    consider_asX.append(0)

        tpcalculus.putSequence(xp, self.filenamex)

        if tpcalculus.exec_():
            print(tpcalculus.getout())

    def saveoutput(self):

        #dir_ = QtWidgets.QFileDialog.getExistingDirectory(None, 'Select a folder:', '.', QtGui.QFileDialog.ShowDirsOnly)

        if self.ui.radioCT.isChecked():
            if self.results.__len__() > self.sequences.__len__():
                print('Multiple calculations...', self.results.__len__(), self.sequences.__len__())

            suffix = '_sCT_'

            strk = "Name" + "\t" + "Edge Number" + "\t" + "Edge connnectivity Index" + "\t" + "Ponderation" + "\t"
            for i in range(1, 16):
                strk = strk + "Momment" + str(i) + "\t"

            strk = strk + '\n'
            name_out = str(self.ui.fileoutput.text()).split('.')
            np.savetxt(name_out[0] + suffix + '.txt', self.results, delimiter='\t', header=strk, fmt='%s')

        else:

            if self.ui.fcm.isChecked():
                suffix = '_fcm_'
                strk = 'Name\tCluster Number\tTotal adjacency index\tPonderation\t'
            else:
                suffix = '_nandy_'
                strk = 'Name\tTotal adjacency index\tPonderation\t'

            if self.ui.linealsort.isChecked():
                suffix = '_lineal_'

            for i in range(1, 16):
                strk = strk + 'Moment ' + str(i) + '\t'

            # print(suffix)
            # print self.results
            # print len(self.results)
            # print self.ponderations
            # print len(self.results[0])

            if self.ponderations.__len__() != 1:
                seedName = str(self.ui.fileoutput.text()).split('.')
                for i in range(self.ponderations.__len__()):
                    name_out = seedName[0] + suffix + str(i) + '.txt'
                    asalvar = []
                    for j in range(self.results.__len__()):

                        if len(self.results[j]) - self.ponderations.__len__() == 1:

                            asalvar.append([self.results[j][0]] + list(self.results[j][i + 1]))
                        elif len(self.results[j]) - self.ponderations.__len__() == 2:

                            asalvar.append([self.results[j][0]] + [self.results[j][1]] + list(self.results[j][i + 2]))

                    np.savetxt(name_out, asalvar, delimiter='\t', header=strk, fmt='%s')
            else:
                name_out = str(self.ui.fileoutput.text()).split('.')

                asalvar = []

                for i in range(self.results.__len__()):
                    if len(self.results[i]) == 2:
                        asalvar.append([self.results[i][0]] + list(self.results[i][1]))
                    elif len(self.results[i]) == 3:
                        asalvar.append([self.results[i][0]] + [self.results[i][1]] + list(self.results[i][2]))

                np.savetxt(name_out[0] + suffix + '.txt', asalvar, delimiter='\t', header=strk, fmt='%s')

        print('Done ...')
        # self.ui.statusbar.showMessage("The calculations were perfomed...")
        # QtGui.QMessageBox.information(self,"Calculation Successfull", "The claculations were performed ...")


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    # app.setStyle('Plastique')
    app.setStyle('cleanlooks')
    app.processEvents()
    myapp = MyInterface()

    # myapp.center()
    myapp.move(QtWidgets.QDesktopWidget().availableGeometry().center() - myapp.frameGeometry().center())

    myapp.show()
    sys.exit(app.exec_())
