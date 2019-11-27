#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
import scipy.cluster.hierarchy as sch
from pylab import savefig
from Bio import SeqIO
import Bio
from itertools import product
import warnings
import re
from  multiprocessing import Pool
from propy import PseudoAAC as PAAC
from propy import Autocorrelation as AC
from propy import CTD as CTD
from propy import QuasiSequenceOrder as QSO
import pseknc_mod as pse

import sys




def dist(X, Y):
    return np.sqrt(np.sum((X - Y) ** 2))

def checknumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def PointDist(point):
    dist = np.zeros([len(point), len(point)])

    for i in range(len(point)):
        for j in range(i, len(point)):
            dist[i][j] = distance.euclidean(point[i], point[j])
            dist[j][i] = dist[i][j]
    return dist


def toworkwith(nps=0):
    if nps == 0:
        seedpattern = Bio.Data.IUPACData.protein_letters
    elif nps == 1:
        seedpattern = 'ACGT'
    elif nps == 2:
        seedpattern = 'ACGU'

    return seedpattern


def CountKmers(listakms, wheretolook):
    # print listakms, wheretolook

    counter = 0
    for kms in listakms:
        counter += re.findall(kms, wheretolook).__len__()
    return counter


def checkPsekNC1(parameters, sequenceX):
    if parameters['Type'] == u'Type I':
        Type = 1
    else:
        Type = 2

    props1 = parameters['Properties']
    lam = parameters['lam']
    kay = parameters['k']
    weight = parameters['w']
    formatt = 'csv'
    geneticMaterial = parameters['geneticMaterial']

    Vector = pse.pseKNC(sequenceX.upper(), props1, Type, kay, lam, weight, formatt, geneticMaterial)

    return Vector


def checkPsekNC2(parameters, sequenceX):
    Vector = pse.simplePseKNC(sequenceX.upper(), int(parameters['k']), 'cvs', parameters['geneticMaterial'])
    return Vector


def checkcomposition(kmerspatrones, musterX, seqforlook):
    result = []
    resultantepatron = []

    # print kmerspatrones
    for x in kmerspatrones:
        resultantepatron.append(replacemuster2(x, musterX))

    for x in resultantepatron:
        result.append(CountKmers(x, seqforlook))

    return result


# PseudoACC module

def calculePseudoAcc(lamdasize, proteinSequence):
    temp = PAAC._GetPseudoAAC(proteinSequence, lamda=lamdasize, weight=0.05)
    return list(temp.values())


# Autocorrelation modules

def calculeNormMoreauBroto(proteinSequence):
    temp = AC.CalculateNormalizedMoreauBrotoAuto(proteinSequence, [AC._ResidueASA], ['ResidueASA'])
    return list(temp.values()[0].values())


def calculeMoran(proteinSequence):
    temp = AC.CalculateMoranAuto(proteinSequence, [AC._ResidueASA], ['ResidueASA'])
    return list(temp.values()[0].values())


def calculeGearyAuto(proteinSequence):
    temp = AC.CalculateGearyAuto(proteinSequence, [AC._ResidueASA], ['ResidueASA'])
    return list(temp.values()[0].values())


def calculateAutoTotal(proteinSequence):
    temp = AC.CalculateAutoTotal(proteinSequence)
    return list(temp.values())


# CTD  module

def calculateC(proteinSequence):
    temp = CTD.CalculateC(proteinSequence)
    return list(temp.values())


def calculateT(proteinSequence):
    temp = CTD.CalculateT(proteinSequence)
    return list(temp.values())


def calculateD(proteinSequence):
    temp = CTD.CalculateD(proteinSequence)
    return list(temp.values())


def calculateCTD(proteinSequence):
    temp = CTD.CalculateCTD(proteinSequence)
    return list(temp.values())


# QuasiSequenceOrder module
def GetSequenceOrderCouplingNumberTotal(proteinSequence, maxlag=30):
    temp = QSO.GetSequenceOrderCouplingNumberTotal(proteinSequence, maxlag=30)
    return list(temp.values())


def GetQuasiSequenceOrder(ProteinSequence, maxlag=30, weight=0.1):
    temp = QSO.GetQuasiSequenceOrder(ProteinSequence, maxlag=30, weight=0.1)
    return list(temp.values())


# calculate spaced kmers
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


def getk_merk(n, NPsort=0):
    seedpattern = toworkwith(NPsort)

    c = []
    c = [''.join(x) for x in product(seedpattern, repeat=n)]

    return c


def checkposX(i, k, l):
    if i.count('X') == k:
        if (i.find('.') != 0) and (i.rfind('.') != k + l - 1):
            return i


class Secuencias:
    def __init__(self):
        self.filename = 'NO_TI2BioP_index'
        self.matrix = []
        self.nmatrix = []
        self.TI2top = []
        self.seqs = []
        self.seqsX = []
        self.seqsMerk = []
        self.PseudoACC = []
        self.PsekNC1 = []
        self.PsekNC2 = []
        self.AutoC = np.array([])
        self.CTD = np.array([])
        self.QuasiSeq = np.array([])
        self.topIndex = False
        self.seqIndex = False


    def readTI2BioP(self, filename):
        self.filename = filename
        f=open(filename, 'r')
        s=f.readline()
        s=f.readline()
        k=s.split()

        aleer = []
        for i in range(len(k)):
            if checknumber(k[i]):
                aleer.append(i)
        f.close()
        print(aleer)
        f=open(filename, 'r')
        s=f.readline()


        #self.TI2top = np.genfromtxt(filename, delimiter='\t', skip_header=1, usecols=tuple([1, 2] + range(4, 19)))
        self.TI2top = np.genfromtxt(filename, delimiter='\t', skip_header=1, usecols=tuple(aleer))

        for i in range(len(self.TI2top)):
            s=f.readline()
            self.seqsX.append(s.split()[0])

        f.close()
        self.topIndex = True


    def readDescript(self, filename):
        self.filename = filename
        f=open(filename, 'r')
        s=f.readline()
        s=f.readline()
        k=s.split(',')
        aleer = []
        for i in range(len(k)):
            if checknumber(k[i]):
                aleer.append(i)
        f.close()
        print(aleer)
        f=open(filename, 'r')
        s=f.readline()
        #self.TI2top = np.genfromtxt(filename, delimiter='\t', skip_header=1, usecols=tuple([1, 2] + range(4, 19)))
        self.TI2top = np.genfromtxt(filename, delimiter=',', skip_header=1, usecols=tuple(aleer))


        for i in range(len(self.TI2top)):
            s=f.readline()
            self.seqsX.append(s.split()[0])

        f.close()
        self.topIndex = True



    def settoCalc(self, n):
        if n == 0:
            self.matrix = self.TI2top
        elif n == 1:
            self.matrix = np.array(self.seqsMerk, dtype=float)
        elif n == 2:
            self.matrix = np.array(self.PseudoACC, dtype=float)
        elif n == 3:
            self.matrix = np.array(self.PsekNC1, dtype=float)
        elif n == 4:
            self.matrix = np.array(self.PsekNC2, dtype=float)
        elif n == 5:
            self.matrix = self.QuasiSeq[0]['Matrix']
        elif n == 6:
            self.matrix = self.QuasiSeq[1]['Matrix']
        elif n == 7:
            self.matrix = self.AutoC[0]['Matrix']
        elif n == 8:
            self.matrix = self.AutoC[1]['Matrix']
        elif n == 9:
            self.matrix = self.AutoC[2]['Matrix']
        elif n == 10:
            self.matrix = self.AutoC[3]['Matrix']
        elif n == 11:
            self.matrix = self.CTD[0]['Matrix']
        elif n == 12:
            self.matrix = self.CTD[1]['Matrix']
        elif n == 13:
            self.matrix = self.CTD[2]['Matrix']
        elif n == 14:
            self.matrix = self.CTD[3]['Matrix']


    def matrixNormalize(self):
        suma = 0.0
        rowsum = []
        sizematrix = len(self.matrix)
        for k in range(sizematrix):
            rowsum.append(self.matrix[k].sum())
            self.matrix[k] = self.matrix[k] / rowsum[k]

    def getSecuenceName(self):
        fichero = open(self.filename, 'r')
        lista = []
        for linea in fichero:
            lista.append(linea.split()[0][1:])
        fichero.close()
        return lista[1:]


    def buildVectorMerk(self, n, NPSort=0):

        Patrones = getk_merk(n, NPSort)

        for i in self.seqsX:
            Spfc = []
            Spfc = [i.seq.count(x) for x in Patrones]
            # print Spfc
            self.seqsMerk.append(Spfc)

    def buildVectorMerk2(self, k, sp, NPSort=0):
        markers = []
        patrones = getk_merk(k, NPSort)
        # print ('I am here.........')
        # print Patrones.__len__()
        muster_template = Gspacedk_mers(k, sp)
        # for x in self.seqsX:
        #     print x.seq
        #     self.seqsMerk.append(checkcomposition(Patrones, muster_template, str(x.seq)))

        pool = Pool()
        resultsx = [pool.apply_async(checkcomposition, (patrones, muster_template, str(x.seq))) for x in self.seqsX]
        [self.seqsMerk.append(r.get()) for r in resultsx]

    def buildPsekNCVectors1(self, parameters):
        pool = Pool()
        resultsx = [pool.apply_async(checkPsekNC1, (parameters, str(x.seq))) for x in self.seqsX]
        [self.PsekNC1.append(r.get()) for r in resultsx]


    def buildPsekNCVectors2(self, parameters):
        pool = Pool()
        resultsx = [pool.apply_async(checkPsekNC2, (parameters, str(x.seq))) for x in self.seqsX]
        [self.PsekNC2.append(r.get()) for r in resultsx]

    def buildvectorsC(self):
        temp = []
        pool2 = Pool()
        resultx = [pool2.apply_async(calculateC, (str(X.seq),)) for X in self.seqsX]
        [temp.append(r.get()) for r in resultx]
        self.CTD = np.append(self.CTD, {'Method': 'C', 'Matrix': np.array(temp)})


    def buildVectorsCTD(self, typeC):
        if typeC == 0:
            temp = []
            pool = Pool()
            resultx = [pool.apply_async(calculateC, (str(x.seq),)) for x in self.seqsX]
            [temp.append(r.get()) for r in resultx]
            self.CTD = np.append(self.CTD, {'Method': 'C', 'Matrix': np.array(temp)})

        elif typeC == 1:
            temp1 = []
            pool2 = Pool()
            resultx = [pool2.apply_async(calculateT, (str(x.seq),)) for x in self.seqsX]
            [temp1.append(r.get()) for r in resultx]
            self.CTD = np.append(self.CTD, {'Method': 'T', 'Matrix': np.array(temp1)})

        elif typeC == 2:
            temp2 = []
            pool3 = Pool()
            resultx = [pool3.apply_async(calculateD, (str(x.seq),)) for x in self.seqsX]
            [temp2.append(r.get()) for r in resultx]
            self.CTD = np.append(self.CTD, {'Method': 'D', 'Matrix': np.array(temp2)})

        elif typeC == 3:
            temp3 = []
            pool4 = Pool()
            resultx = [pool4.apply_async(calculateCTD, (str(x.seq),)) for x in self.seqsX]
            [temp3.append(r.get()) for r in resultx]
            self.CTD = np.append(self.CTD, {'Method': 'CTD', 'Matrix': np.array(temp3)})

    def buildVectorsAutoC(self, typeA):
        if typeA == 0:
            temp1 = []
            pool1 = Pool()
            resultx = [pool1.apply_async(calculeNormMoreauBroto, (str(x.seq),)) for x in self.seqsX]
            [temp1.append(r.get()) for r in resultx]
            self.AutoC = np.append(self.AutoC, {'Method': 'NormMoreauBroto', 'Matrix': np.array(temp1)})
        elif typeA == 1:
            temp2 = []
            pool2 = Pool()
            resultx = [pool2.apply_async(calculeMoran, (str(x.seq),)) for x in self.seqsX]
            [temp2.append(r.get()) for r in resultx]
            self.AutoC = np.append(self.AutoC, {'Method': 'Moran', 'Matrix': np.array(temp2)})
        elif typeA == 2:
            temp3 = []
            pool3 = Pool()
            resultx = [pool3.apply_async(calculeGearyAuto, (str(x.seq),)) for x in self.seqsX]
            [temp3.append(r.get()) for r in resultx]
            self.AutoC = np.append(self.AutoC, {'Method': 'Geary', 'Matrix': np.array(temp3)})
        elif typeA == 3:
            temp4 = []
            pool4 = Pool()
            resultx = [pool4.apply_async(calculateAutoTotal, (str(x.seq),)) for x in self.seqsX]
            [temp4.append(r.get()) for r in resultx]
            self.AutoC = np.append(self.AutoC, {'Method': 'AutoTotal', 'Matrix': np.array(temp4)})

    def buildVectorsQuasi(self, typeQ):
        if typeQ['t'] == 0:
            tempx = []
            pool = Pool()
            resultx = [pool.apply_async(GetSequenceOrderCouplingNumberTotal, ((str(x.seq), typeQ['Maxlag']))) for x in
                       self.seqsX]
            [tempx.append(r.get()) for r in resultx]
            self.QuasiSeq = np.append(self.QuasiSeq, {'Method': 'Quasi_Sequence_Coupling', 'Matrix': np.array(tempx)})
        elif typeQ['t'] == 1:
            temp1 = []
            pool1 = Pool()
            resultx = [
                pool1.apply_async(GetQuasiSequenceOrder, (str(x.seq), float(typeQ['Maxlag']), float(typeQ['Weight'])))
                for x in self.seqsX]
            [temp1.append(r.get()) for r in resultx]
            self.QuasiSeq = np.append(self.QuasiSeq, {'Method': 'Quasi_Sequence_Order', 'Matrix': np.array(temp1)})


    def savetolog(self, LogFile, select, spaced=[]):
        LF = open(LogFile, "a")
        np.set_printoptions(precision=3)
        np.set_printoptions(suppress=True)

        if select == 0:
            if len(spaced) != 0:
                LF.write('\nKmer-vectors: k= ' + str(spaced[0])+' spaced= '+str(spaced[1]) + '\n')
            else:
                LF.write('\nKmer-vectors:\n')
            for x in range(len(self.seqsMerk)):
                LF.write(str(self.seqsX[x].id) + ' --> ' + str(self.seqsMerk[x]) + '\n')

        elif select == 1:
            LF.write('\nPseudoACC-vectors:\n')
            for x in range(len(self.seqsX)):
                LF.write(str(self.seqsX[x].id) + '--> ' + str(self.PseudoACC) + '\n')
        elif select == 2:
            LF.write('\nPsekNC-vectors:\n')
            for x in range(len(self.seqsX)):
                LF.write(str(self.seqsX[x].id) + '--> ' + str(self.PsekNC1[x]) + '\n')

        elif select == 3:
            LF.write('\nSimple PsekNC-vectors:\n')
            for x in range(len(self.seqsX)):
                LF.write(str(self.seqsX[x].id) + '--> ' + str(self.PsekNC2[x]) + '\n')
        elif select == 4:
            LF.write('\n CTD-Vectors')
            for x in self.CTD:
                LF.write('\n' + str(x['Method']) + '-vectors:\n')
                for y in range(len(self.seqsX)):
                    LF.write(self.seqsX[y].id + '--> ' + str(x['Matrix'][y]) + '\n')
        elif select == 5:
            LF.write('\n Autocorrection-vectors:')
            for x in self.AutoC:
                LF.write('\n' + str(x['Method']) + '-vectors:\n')
                for y in range(len(self.seqsX)):
                    LF.write(self.seqsX[y].id + '--> ' + str(x['Matrix'][y]) + '\n')
        elif select == 6:
            LF.write('\n Quasi Sequence Order-vectors:')
            for x in self.QuasiSeq:
                LF.write('\n' + str(x['Method']) + '-vectors:\n')
                for y in range(len(self.seqsX)):
                    LF.write(self.seqsX[y].id + '--> ' + str(x['Matrix'][y]) + '\n')

        LF.close()


    def buildVectorPseudoACC(self, lamdasize):
        pool = Pool()
        resultsx = [pool.apply_async(calculePseudoAcc, (lamdasize, str(x.seq))) for x in self.seqsX]
        [self.PseudoACC.append(r.get()) for r in resultsx]


    def ShowDendrogram(self):
        d = sch.distance.pdist(self.matrix)
        print (d)
        Z = sch.linkage(d, method='complete')
        P = sch.dendrogram(Z, orientation='right', labels=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'])
        savefig('testing.pdf')
        plt.show()


    def ShowDendrogramChoice(self, DistSelected=0, order=0, prob=0):
        if DistSelected == 0:
            d = sch.distance.pdist(self.matrix)
        elif DistSelected == 1:
            d = self.BrasCurtisDistance()
        elif DistSelected == 2:
            d = self.HammingDistance()
        elif DistSelected == 3:
            d = self.CanberraDistance()
        elif DistSelected == 4:
            d = self.ChebyshevDistance()
        elif DistSelected == 5:
            d = self.MinkowskiDistance(p=prob)
        elif DistSelected == 6:
            d = self.WminkowskiDistance(order)
        elif DistSelected == 7:
            d = self.Jensen_ShannonDistance()

        Z = sch.linkage(d, method='complete')
        P = sch.dendrogram(Z, orientation='right', labels=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'])
        savefig('testing.pdf')
        plt.show()


    def EuclidianDistance(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.euclidean(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]

        return DistMat


    def WminkowskiDistance(self, p):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        vector = np.ones(n)
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.wminkowski(self.matrix[i], self.matrix[j], p, vector)
                DistMat[j, i] = DistMat[i, j]

        return DistMat


    def MinkowskiDistance(self, p):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        # vector = np.ones(15) # hay versiones que necesitan el vector unitario
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.minkowski(self.matrix[i], self.matrix[j], p)
                DistMat[j, i] = DistMat[i, j]

        return DistMat


    def Kulsinski(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.kulsinski(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]
        return DistMat


    def HammingDistance(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.hamming(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]
        return DistMat


    def BrasCurtisDistance(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        for i in range(n):
            for j in range(i + i, n):
                DistMat[i, j] = distance.braycurtis(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]
        return DistMat


    def CanberraDistance(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.canberra(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]
        return DistMat


    def ChebyshevDistance(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        vector = np.ones(15)
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = distance.chebyshev(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]

        return DistMat


    def jsd(self, x, y):
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        x = np.array(x)
        y = np.array(y)
        d1 = x * np.log2(2 * x / (x + y))
        d2 = y * np.log2(2 * y / (x + y))
        d1[np.isnan(d1)] = 0
        d2[np.isnan(d2)] = 0
        return 0.5 * np.sum(d1 + d2)


    def Jensen_ShannonDistance(self):
        n = self.matrix.__len__()
        DistMat = np.zeros([n, n])
        vector = np.ones(15)
        for i in range(n):
            for j in range(i + 1, n):
                DistMat[i, j] = self.jsd(self.matrix[i], self.matrix[j])
                DistMat[j, i] = DistMat[i, j]

        return DistMat


    def SeqListVector(self, filename, NPSort=0):
        aminolist = {}
        allprob = []
        self.seqsX = list(SeqIO.parse(filename, 'fasta'))
        print(self.seqsX)


        for x in self.seqsX:
            x.seq = x.seq.upper()

        seekpattern = toworkwith(NPSort)

        # AminoList = Bio.Data.IUPACData.protein_letters
        for sequence in self.seqsX:
            aminolist = []
            SeqSize = float(len(sequence))

            for AminoAcid in seekpattern:
                aminolist.append(sequence.seq.count(AminoAcid) / SeqSize)
            allprob.append(aminolist)

        return allprob


    def writeMega(self, distMatrix, filename_out, matrixName):
        fichero = open(filename_out, 'w')
        print ('writing .... ' + filename_out)
        slist = '[' + '     '
        lprintd = '['
        fichero.write('#mega \n')

        fichero.write('!Title: ' + matrixName + '\t\n')
        fichero.write(';\t\n')
        fichero.write('!Format DataType=distance; \n')

        if self.seqsX:
            for element in self.seqsX:
                fichero.write('#' + element.id + '\n')
                slist = slist + element.id + '       '
            fichero.write('\n')
            fichero.write(slist + ']' + '\n')
        for i in range(len(distMatrix)):
            lprintd = '[' + self.seqsX[i].id + ']' + '  '
            for j in range(i):
                lprintd = lprintd + str(distMatrix[i][j]) + '\t'

            fichero.write(lprintd + '\n')
        fichero.close()


if __name__ == '__main__':
    secuencia = Secuencias('seqs_0T_out.txt')
    secuencia.settoCalc(0)
    secuencia.ShowDendrogramChoice(DistSelected=5, prob=0.5)
    secuencia.seqs = secuencia.SeqListVector('seqs_0T.txt')
    secuencia.settoCalc(1)
    secuencia.buildVectorMerk(4)
    print(secuenciencia.seqsMerk)
    d = np.zeros((10, 15))
    secuencia.writeMega(d, 'testing_xx.txt', 'distance_matrix_0T')
    secuencia.ShowDendrogramChoice(DistSelected=7)


