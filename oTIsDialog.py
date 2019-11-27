__author__ = 'reymolina'

import sys
import re
import Bio
import numpy as np

import FastaProc

from PyQt5 import QtGui, QtCore, QtWidgets
from  multiprocessing import Pool
from itertools import product
from propy import PseudoAAC as PAAC
from propy import Autocorrelation as AC
from propy import CTD as CTD
from propy import QuasiSequenceOrder as QSO
import pseknc_mod as pse

from TIs4AAQt import *


# Relative composition

def relativecomposition(name, pattern, secuenciax ):
    aminolist = {}
    aminolist['name'] = name
    seqsize = float(secuenciax.__len__())
    for aminoacid in pattern:
        aminolist[aminoacid] = secuenciax.count(aminoacid)/ seqsize

    return  aminolist


# PseudoACC module

def calculePseudoAcc(lamdasize, proteinSequence, name):
    # return a dictionary with PseudoAAC
    temp = PAAC._GetPseudoAAC(proteinSequence, lamda=lamdasize, weight=0.05)
    temp ['name'] = name
    return temp


# Autocorrelation modules

def calculeNormMoreauBroto(proteinSequence, name):
    temp = AC.CalculateNormalizedMoreauBrotoAuto(proteinSequence, [AC._ResidueASA], ['ResidueASA'])
    temp['name'] = name
    return temp

def calculeMoran(proteinSequence, name):
    temp = AC.CalculateMoranAuto(proteinSequence, [AC._ResidueASA], ['ResidueASA'])
    temp['name'] = name
    return temp

def calculeGearyAuto(proteinSequence, name):
    temp = AC.CalculateGearyAuto(proteinSequence, [AC._ResidueASA], ['ResidueASA'])
    temp['name'] = name
    return temp

def calculateAutoTotal(proteinSequence, name):
    temp = AC.CalculateAutoTotal(proteinSequence)
    temp ['name'] = name
    return temp


# CTD  module

def calculateC(proteinSequence, name):
    temp = CTD.CalculateC(proteinSequence)
    temp['name'] = name
    return temp


def calculateT(proteinSequence, name):
    temp = CTD.CalculateT(proteinSequence)
    temp['name'] = name
    return temp


def calculateD(proteinSequence, name):
    temp = CTD.CalculateD(proteinSequence)
    temp['name'] = name
    return temp


def calculateCTD(proteinSequence, name):
    temp = CTD.CalculateCTD(proteinSequence)
    temp['name'] = name
    return temp


# QuasiSequenceOrder module

def GetSequenceOrderCouplingNumberTotal(proteinSequence, name, maxlag=30):
    temp = QSO.GetSequenceOrderCouplingNumberTotal(proteinSequence, maxlag=30)
    temp['name'] = name
    return temp


def GetQuasiSequenceOrder(ProteinSequence, name, maxlag=30, weight=0.1):
    temp = QSO.GetQuasiSequenceOrder(ProteinSequence, maxlag=30, weight=0.1)
    temp['name'] = name
    return temp


# K-mers without spaces

def getk_merk(n, seedpattern):
    c = []
    c = [''.join(x) for x in product(seedpattern, repeat=n)]
    return c


def GetKmersSeq(Protein, name, patrones):
    kmers = {}
    kmers['name'] = name
    for x in patrones:
        kmers[x] = Protein.count(x)
    return kmers


# Spaced Kmers
# calculate spaced kmers


def checkposX(i, k, l):
    if i.count('X') == k:
        if (i.find('.') != 0) and (i.rfind('.') != k + l - 1):
            return i


def Gspacedk_mers(k, l):
    result = []
    istring = k * 'X' + l * '.'
    c = [''.join(x) for x in product(istring, repeat=k + l)]
    [result.append(checkposX(x, k, l)) for x in c]
    result = list(set(result))
    del (result[result.index(None)])

    return result


def replacemuster2(patroneselement, musterlist):
    resultsf = [replaceMuster(patroneselement, x) for x in musterlist]
    return resultsf


def replaceMuster(patronelement, muster):
    cad = ''
    counter = 0
    for j in range(len(muster)):
        if muster[j] == 'X':
            cad = cad + patronelement[counter]
            counter += 1
        else:
            cad = cad + '.'

    return cad


def CountKmers(listakms, wheretolook):
    # print listakms, wheretolook

    counter = 0
    for kms in listakms:
        counter += re.findall(kms, wheretolook).__len__()
    return counter

def checkcomposition(kmerspatrones, musterX, seqforlook, name):
    result = {}
    result['name'] = name
    resultantepatron = []

    # print kmerspatrones
    for x in kmerspatrones:
        resultantepatron.append(replacemuster2(x, musterX))

    for x in resultantepatron:
        result[x[0].replace('.','0')] = CountKmers(x, seqforlook)

    return result

def saveresultx(filename, datax):
    header = datax[0].keys()
    header.remove('name')
    header.insert(0,'name')
    f = open(filename, 'w')
    strtosave = ''
    for x in header:
        strtosave += ',' + x

    f.writelines(strtosave + '\n')
    for x in datax:
        strtosave = ''
        for keyX in header:
            strtosave += ',' + str(x[keyX])

        f.writelines(strtosave + '\n')


    f.close()


class DialogTIs(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self,parent)
        self.ui=Ui_TIsforAA()
        self.ui.setupUi(self)
        self.calcule = False
        self.behavior = 'Maybe next time...'
        self.ui.Select_all.clicked.connect(self.SelectAll)
        self.ui.tocalc.clicked.connect(self.calculate_all)
        self.seedpattern = Bio.Data.IUPACData.protein_letters

    def SelectAll(self):
        self.ui.checkBox.setChecked(True)
        self.ui.checkBox_2.setChecked(True)

        self.ui.pseudoCalc.setChecked(True)

        self.ui.checkBox_9.setChecked(True)
        self.ui.checkBox_10.setChecked(True)
        self.ui.checkBox_11.setChecked(True)
        self.ui.checkBox_12.setChecked(True)

        self.ui.checkBox_13.setChecked(True)
        self.ui.checkBox_14.setChecked(True)
        self.ui.checkBox_15.setChecked(True)
        self.ui.checkBox_16.setChecked(True)

    def calculate_all(self):

        if self.calcule:
            self.ui.progressBar.setRange(0, 0)
            self.ui.progressBar.setValue(0)

            if self.ui.pseudoCalc.isChecked():
                print ('Calculation of Pseudo Amino Acid composition ...')

                lamdax = int(self.ui.spinBox_4.value())
                filenamex = self.filenamex + '_PseudoACC_Lamda_' + str(lamdax) + '.csv'
                print('lamdax = ', lamdax)
                self.PseudoAcc = []
                pool = Pool()
                resultsx = [pool.apply_async(calculePseudoAcc, (lamdax, x[1], x[0])) for x in self.sequences]
                [self.PseudoAcc.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.PseudoAcc)
                self.PseudoAcc = []

                # for x in self.sequences:
                #     print ('{} -> secuencia:{}'.format(x[0], x[1]))

            filenamex = self.filenamex + '_relative_composition.csv'
            self.composition = []
            print ('Relative amino acid composition ...')
            pool = Pool()
            resultsx = [pool.apply_async(relativecomposition, (x[0], self.seedpattern, x[1])) for x in self.sequences]
            [self.composition.append(r.get()) for r in resultsx]
            saveresultx(filenamex, self.composition)
            self.composition=[]

            #Kmers without space

            self.Kmers = []
            print ('Calculation of  Kmers ....' )
            sizeKmers = int(self.ui.spinBox_2.value())
            filenamex = self.filenamex + '_' + str(sizeKmers) + '_kmers.csv'
            sizeSpace = int(self.ui.spinBox.value())
            patrones = getk_merk(sizeKmers, self.seedpattern)
            pool = Pool()
            resultsx = [pool.apply_async(GetKmersSeq, (x[1], x[0], patrones )) for x in self.sequences]
            [self.Kmers.append(r.get()) for r in resultsx]
            saveresultx(filenamex, self.Kmers)
            self.Kmers = []

            if sizeSpace != 0:
                print('Calculation of  Spaced Kmers ....')
                muster_template = Gspacedk_mers(sizeKmers, sizeSpace)
                self.spacedKmers = []
                pool = Pool()
                resultsx = [pool.apply_async(checkcomposition, (patrones, muster_template, x[1], x[0])) for x in self.sequences]
                [self.spacedKmers.append(r.get()) for r in resultsx]
                filenamex = self.filenamex + '_' + str(sizeKmers)+ '_kmers_space_' + str(sizeSpace) + '.csv'
                saveresultx(filenamex,self.spacedKmers)
                self.spacedKmers = []

            if  self.ui.checkBox_12.isChecked():
                print ('calulation of total autocorrelation ...')
                filenamex = self.filenamex + '_Total_autocorrelation.csv'
                self.Autototal = []
                pool = Pool()
                resultsx = [pool.apply_async(calculateAutoTotal, (x[1], x[0])) for x in self.sequences]
                [self.Autototal.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.Autototal)
                self.Autototal = []

            if  self.ui.checkBox_11.isChecked():
                print ('calulation of Geary autocorrelation ...')
                filenamex = self.filenamex + '_Geary_autocorrelation.csv'
                self.GearyAuto = []
                pool = Pool()
                resultsx = [pool.apply_async(calculeGearyAuto, (x[1], x[0])) for x in self.sequences]
                [self.GearyAuto.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.GearyAuto)
                self.GearyAuto = []

            if  self.ui.checkBox_10.isChecked():
                print ('calulation of Moran autocorrelation ...')
                filenamex = self.filenamex + '_Moran_autocorrelation.csv'
                self.MoranAuto = []
                pool = Pool()
                resultsx = [pool.apply_async(calculeMoran, (x[1], x[0])) for x in self.sequences]
                [self.MoranAuto.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.MoranAuto)
                self.MoranAuto = []


            if  self.ui.checkBox_9.isChecked():
                print ('calulation of Moreau Broto autocorrelation ...')
                filenamex = self.filenamex + '_MoreauBrotoAuto.csv'
                self.MoreauBrotoAuto = []
                pool = Pool()
                resultsx = [pool.apply_async(calculeNormMoreauBroto, (x[1], x[0])) for x in self.sequences]
                [self.MoreauBrotoAuto.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.MoreauBrotoAuto)
                self.MoreauBrotoAuto = []


            if  self.ui.checkBox_13.isChecked():
                print ('calulation of C from CTD ...')
                filenamex = self.filenamex + '_CTDC.csv'
                self.CTDC = []
                pool = Pool()
                resultsx = [pool.apply_async(calculateC, (x[1], x[0])) for x in self.sequences]
                [self.CTDC.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.CTDC)
                self.CTDC = []

            if  self.ui.checkBox_14.isChecked():
                filenamex = self.filenamex + '_CTDT.csv'
                print ('calulation of T from CTD ...')
                self.CTDT = []
                pool = Pool()
                resultsx = [pool.apply_async(calculateT, (x[1], x[0])) for x in self.sequences]
                [self.CTDT.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.CTDT)
                self.CTDT = []

            if  self.ui.checkBox_15.isChecked():
                filenamex = self.filenamex + '_CTDD.csv'
                print ('calulation of D from CTD ...')
                self.CTDD = []
                pool = Pool()
                resultsx = [pool.apply_async(calculateD, (x[1], x[0])) for x in self.sequences]
                [self.CTDD.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.CTDD)
                self.CTDD = []

            if  self.ui.checkBox_16.isChecked():
                filenamex = self.filenamex + '_CTD.csv'
                print ('calulation of CTD  ...')
                self.CTD = []
                pool = Pool()
                resultsx = [pool.apply_async(calculateCTD, (x[1], x[0])) for x in self.sequences]
                [self.CTD.append(r.get()) for r in resultsx]
                saveresultx(filenamex, self.CTD)
                self.CTD = []

            if self.ui.checkBox.isChecked() or self.ui.checkBox_2.isChecked():
                maxlagx = int(self.ui.spinBox_3.value())
                weightx = float(self.ui.lineEdit.text())

                if self.ui.checkBox.isChecked():
                    print ('calulation of Sequence Order Coupling ...')
                    filenamex = self.filenamex + '_SQOCoupling.csv'
                    self.Coupling = []
                    pool = Pool()
                    resultx = [pool.apply_async(GetSequenceOrderCouplingNumberTotal, (x[1], x[0], maxlagx)) for x in self.sequences]
                    [self.Coupling.append(r.get()) for r in resultsx]
                    saveresultx(filenamex, self.Coupling)
                    self.Coupling = []

                if self.ui.checkBox_2.isChecked():
                    print ('calulation of Quasi Sequence Order ...')
                    filenamex = self.filenamex + '_QuasiSQO.csv'
                    self.Quasi = []
                    pool = Pool()
                    resultx = [pool.apply_async(GetQuasiSequenceOrder, (x[1], x[0], maxlagx, weightx)) for x in self.sequences]
                    [self.Quasi.append(r.get()) for r in resultsx]
                    saveresultx(filenamex, self.Quasi)
                    self.Quasi = []



            self.behavior = 'Calculation were performed ....'





            print ('... Done.')
            self.ui.progressBar.setRange(0, 1)
            self.close()
        else:
            print (' nothing to do ....')

    def getout(self):
        return self.behavior

    def putSequence(self, sequences, filenamex):

        self.sequences = sequences
        self.filenamex = filenamex
        if self.sequences.__len__() != 0:
            self.calcule = True
        else:
            self.calcule = False








if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = DialogTIs()
    myapp.Addponderations()
    if myapp.exec_():
        values = myapp.getvalues()
        print(values)


    #myapp.show()

    #sys.exit(app.exec_())

