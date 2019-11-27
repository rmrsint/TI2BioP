import sys
import datapond
import numpy as np
from Bio.SeqUtils import seq1
from PyQt5 import QtGui, QtCore, QtWidgets
from PIL import Image, ImageDraw
from multiprocessing import Pool

orientaseq = ('G', 'C', 'A', ['T', 'U'])
orientAmino = (
    ['A', 'L', 'V', 'I', 'P', 'F', 'M', 'W'], ['G', 'S', 'T', 'N', 'Q', 'Y', 'C'], ['D', 'E', 'B', 'Z'],
    ['K', 'R', 'H'])
orienta_seq_lineal = ('', '', '', ['C', 'G', 'A', 'T', 'U'])
orienta_amino_lineal = ('', '', '',
                        ['A', 'L', 'V', 'I', 'P', 'F', 'M', 'W', 'G', 'S', 'T', 'N', 'Q', 'Y', 'C', 'D', 'E', 'B', 'Z',
                         'K', 'R', 'H'])

def totaladj_index(mcontactos):
    counter = 0
    for i in range(len(mcontactos)):
        for j in range(len(mcontactos)):
            if j > i and mcontactos[i][j] == 1:
                counter += 1
    return 2 * counter

class MyWidget(QtWidgets.QWidget):
    def __init__(self, distribution, bonds,  title=None,  savefilename='mydrawingx.jpg'):
        super(MyWidget, self).__init__()
        self.distribution = distribution
        self.bonds = bonds
        self.title = title
        if '>' in savefilename:
            self.fname = savefilename[1:]
        else:
            self.fname = savefilename
        self.initUI()
        self.test = self
        self.label = QtWidgets.QLabel(self)

    def initUI(self):
        self.setGeometry(10, 10, 800, 900)
        self.setWindowTitle(self.title)
        self.show()


    def paintEvent(self, e):
        qp = QtGui.QPainter()
        self.drawNandy(qp)


    def drawNandy(self, qp):
        nandydata = self.distribution
        bonds = self.bonds
        # print nandydata
        n = nandydata.__len__()
        listado = {}
        w = 800
        h = 700
        sizex = 50
        sizey = 50
        xlower = 1000
        ylower = 1000
        xmax = -3200
        ymax = -3200
        pixmap = QtGui.QPixmap(w, h)
        pixmap.fill(QtCore.Qt.white)
        qp.begin(pixmap)
        qp.setBrush(QtGui.QColor(255, 255, 255))  # white was 204, 204, 0
        qp.drawRect(0, 0, w, h)
        qp.setBrush(QtGui.QColor(221, 34, 0)) # was 221, 34, 0

        for x in nandydata:
            listado[nandydata[x][1]] = x

        for x in bonds:

            if listado[x[0]][0] > xmax:
                xmax = listado[x[0]][0]
            elif listado[x[0]][0] < xlower:
                xlower = listado[x[0]][0]

            if listado[x[1]][0] > xmax:
                xmax = listado[x[1]][0]
            elif listado[x[1]][0] < xlower:
                xlower = listado[x[1]][0]

            if listado[x[0]][1] > ymax:
                ymax = listado[x[0]][1]
            elif listado[x[0]][1] < ylower:
                ylower = listado[x[0]][1]

            if listado[x[1]][1] > ymax:
                ymax = listado[x[1]][1]
            elif listado[x[1]][1] < ylower:
                ylower = listado[x[1]][1]

        offsetx = w // 2 - (xlower + xmax) // 2
        offsety = h // 2 - (ylower + ymax) // 2
        print('lowerx = {}, xmaxm = {}, ylower = {} , ymax = {}, '.format(xlower, xmax, ylower, ymax))
        print('offsetx = {}, offsety = {}'.format(offsetx, offsety))

        factor = offsetx / offsety
        if factor > 1.0:
            factor -= 1
        elif factor < 0.5:
            factor = 0.5

        if (xlower > 0) and (xmax < w):
            if (ylower > 0) and (ymax < h):
                factor = 1.0

        print('factor = {}'.format(factor))
        pen = QtGui.QPen()
        qp.setBrush(QtGui.QColor(220, 34, 0))

        pen.setColor(QtGui.QColor(20, 0, 255))
        pen.setWidth(3)
        qp.setPen(pen)

        makecenter = False
        for x in bonds:
            x0 = listado[x[0]][0]
            y0 = listado[x[0]][1]
            x1 = listado[x[1]][0]
            y1 = listado[x[1]][1]

            if not makecenter:
                qp.drawLine((x0 + offsetx) * factor, 0, (x0 + offsetx) * factor, h)
                qp.drawLine(0, (y0 + offsety) * factor, w, (y0 + offsety)*factor)
                makecenter = True

            pen.setColor(QtGui.QColor(0, 0, 0))
            pen.setWidth(3)
            qp.setPen(pen)
            qp.drawLine((x0 + offsetx) * factor, (y0 + offsety) * factor, (x1 + offsetx) * factor, (y1 + offsety) * factor)
            pen.setWidth(5)
            qp.setPen(pen)
            qp.setBrush(QtGui.QColor(0, 0, 0))
            rectangle = QtCore.QRect((x0 + offsetx) * factor, (y0 + offsety) * factor, 1, 1)
            qp.drawEllipse(rectangle)
            rectangle = QtCore.QRect((x1 + offsetx) * factor, (y1 + offsety) * factor, 1, 1)
            qp.drawEllipse(rectangle)

        qp.setBrush(QtGui.QColor(21, 34, 0))
        rectangle = QtCore.QRect(w / 2 + offsetx, h / 2 + offsety, 1, 1)
        qp.drawEllipse(rectangle)
        qp.end()

        self.label.setPixmap(pixmap)
        self.label.show()
        self.label.pixmap().save(self.fname)


class MacroMol():
    def __init__(self, seqraw, orientation, ponderations, seqsort=0, show=False, calc=True):
        self.sequence = ''
        self.names = seqraw[0]
        self.seqarray = []
        self.topindexs = []
        self.macrotyp = 0
        self.orientation = orientation
        self.ponderation = ponderations
        self.distribution = {}
        self.bonds = {}
        self.result = []
        self.result.append(self.names)
        self.all_in_one(seqraw[1], seqsort, show=show, perfcalc=calc)

    def addsequence(self, s, seqtyp):
        if seqtyp == 0:
            self.macrotyp = 0
            self.sequence = s.upper()
        elif seqtyp == 1:
            self.macrotyp = 0
            self.sequence = seq1(s).upper()
        elif seqtyp == 2:
            self.macrotyp = 1
            self.sequence = s.upper()
        elif seqtyp == 3:
            self.macrotyp = 2
            self.sequence = s.upper()

        self.n = self.sequence.__len__()

    def nandydatagraph(self):
        up_seq = self.orientation[0]
        down_seq = self.orientation[1]
        left_seq = self.orientation[2]
        right_seq = self.orientation[3]
        seq_scan = self.sequence[1:]
        posx = 150
        posy = 150
        counter = 1
        cluster = []
        bonds = []
        step = 25
        clusternumber = 0
        distribution = {(posx, posy): [self.sequence.upper()[0] + '1', 0]}
        cluster.append((posx, posy))
        previus_pos = 0

        for element in seq_scan:
            counter += 1

            if element in up_seq:
                posy = posy - step
            elif element in down_seq:
                posy = posy + step
            elif element in right_seq:
                posx = posx + step
            elif element in left_seq:
                posx = posx - step

            if (posx, posy) in distribution:
                distribution[(posx, posy)][0] = distribution[(posx, posy)][0] + ';' + element + str(counter)
                current_pos = int(distribution[(posx, posy)][1])
            else:
                distribution[(posx, posy)] = [element + str(counter), 0]
                clusternumber += 1
                current_pos = clusternumber
                distribution[(posx, posy)][1] = clusternumber

            cluster.append((posx, posy))

            if current_pos != previus_pos:

                if current_pos < previus_pos:
                    if not (current_pos, previus_pos) in bonds:
                        bonds.append((current_pos, previus_pos))
                else:
                    if not (previus_pos, current_pos) in bonds:
                        bonds.append((previus_pos, current_pos))

            previus_pos = current_pos

        self.distribution = distribution
        self.bonds = bonds

    def getdata(self):
        return self.distribution, self.bonds


    def buildgraph(self):
        nandydata = self.distribution
        bonds = self.bonds

        n = nandydata.__len__()
        listado = {}
        # 65,100,137,0

        #im = Image.new("RGBA", (300, 300), color=(255, 255, 255, 0))
        im = Image.new("RGBA", (800, 600), color=(24, 231, 24, 0))

        Draw = ImageDraw.Draw(im)

        for x in nandydata:
            listado[nandydata[x][1]] = x
        for x in bonds:
            print(x)
            Draw.line([listado[x[0]], listado[x[1]]], fill=(65, 100, 137), width=2)

        return im

    def showNandy(self):
        presentar = MyWidget(self.distribution, self.bonds, title='Nandy representation' + '  ' + self.names, savefilename=self.names)

        #sys.exit(app.exec_())


    def setPond(self, ponddata):
        pondcluster = np.zeros(self.distribution.__len__())
        if self.macrotyp == 0:
            id = datapond.AASEQ2
        else:
            id = datapond.NAcids
        qcluster = self.distribution.values()
        for x in qcluster:
            for z in x[0].split(';'):
                index = id.find(z[0])
                pondcluster[x[1]]= pondcluster[x[1]] + ponddata[index]
        return pondcluster

    def buildconnectivity(self):
        self.connectivity=np.zeros([self.distribution.__len__(), self.distribution.__len__()])
        for i in self.bonds:
            x, y  = i
            self.connectivity[x][y] = 1
            self.connectivity[y][x] = 1
        for j in range(self.connectivity.__len__()):
            self.connectivity[j][j] = 1



    def calc(self):
        self.buildconnectivity()
        totaladj = totaladj_index(self.connectivity)
        for k in self.ponderation:
            result = []
            pondname, data = datapond.getPond(k, self.macrotyp)
            result.append(totaladj)
            result.append(pondname)
            ponderar = self.setPond(data)

            matrixcalc = np.zeros([len(self.connectivity), len(self.connectivity)])

            for i in range(len(self.connectivity)):
                for j in range(len(self.connectivity)):
                    if j >= i:
                        if self.connectivity[i][j] == 1:
                            matrixcalc[i][j] = (ponderar[i] + ponderar[j] ) / 2
                            matrixcalc[j][i] = matrixcalc[i][j]
            sumcol = np.zeros(len(matrixcalc))

            for i in range(len(matrixcalc)):
                sumcol[i] = matrixcalc[i][:].sum()
                if sumcol[i]== 0:
                    sumcol[i] = 1.0


            matrixcalc = matrixcalc / sumcol
            auxmatrix = matrixcalc.copy()
            result.append(round(matrixcalc.trace(), 4))

            for i in range(14):
                auxmatrix = auxmatrix.dot(matrixcalc)
                result.append(round(auxmatrix.trace(), 4))

            self.result.append(result)


    def all_in_one(self, sequence, ksort, show=False, perfcalc=False):
        self.addsequence(sequence, ksort)
        # print self.n

        self.nandydatagraph()
        if show:
            self.showNandy()

        if perfcalc:
            self.calc()

def calcdescriptors(seqdata, orientation, ponderations, seqsort):
    mol = MacroMol(seqdata, orientation, ponderations, seqsort, calc=True)

    return mol.result

def process_all(sequences, orientation, ponderations, seqsort):
    pool = Pool()
    temp = []
    restx = [pool.apply_async(calcdescriptors, (x, orientation, ponderations, seqsort, ), )
    for x in sequences]

    [temp.append(r.get()) for r in restx]

    return temp

def readfilefasta(filename):
    f = open(filename, 'r')
    sequences = f.readlines()
    depurated = []

    #remove empty lines
    #print len(sequences)
    for x in sequences:
        if x == '\n':
            sequences.remove('\n')

    nspliter = len(sequences) / 2
    for i in range(nspliter):
        depurated.append([sequences[2*i], sequences[2*i +1]])

    return depurated

if __name__ == "__main__":

    macromolecule = MacroMol(['> xx1',
                              'ccgcgcgtgtgtgtgtgtgtgcgctcgctcgctgagcgatcgcatacgactcagcatcagcagcgcatacgcatcagcatcagcatcagcagcatcagcatcagactacgcgcgatatactcagtactgcagtcagtagtcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctatatatatgctacgtacgcatgactcagactaccg'],
                             orientaseq, [2], 2, False, True)
    print (macromolecule.result)

    prot = 'MKHFSKLCFLLSTFAVSIAPVTWAHEGATHQHANVSKLTDAYTYANYDQVKATHVYLDLNVDFDKKSLSG' \
           'FAELSLDWFTDNKAPLILDTRDLVIHRVMAKNSQGQWVKVNYDLAKRDDVLGSKLTINTPLNAKKVRVYY' \
           'NSTEKATGLQWLSAEQTAGKEKPFLFSQNQAIHARSWIPIQDTPSVRVTYTARITTDKDLLAVMSANNEP' \
           'GTERDGDYFFSMPQAIPPYLIAIGVGDLEFKAMSHQTGIYAESYILDAAVAEFDDTQAMIDKAEQMYGKY' \
           'RWGRYDLLMLPPSFPFGGMENPRLSFITPTVVAGDKSLVNLIAHELAHSWSGNLVTNESWRDLWLNEGFT' \
           'SYVENRIMEAVFGTDRAVMEQALGAQDLNAEILELDASDTQLYIDLKGRDPDDAFSGVPYVKGQLFLMYL' \
           'EEKFGRERFDAFVLEYFDSHAFQSLGTDNFVKYLKANLTDKYPNIVSDNEINEWIFKAGLPSYAPQPTSN' \
           'AFKVIDKQINQLVTDELTLEQLPTAQWTLHEWLHFINNLPVDLDHQRMVNLDKAFDLTNSSNAEIAHAWY' \
           'LLSVRADYKEVYPAMAKYLKSIGRRKLIVPLYKELAKNAESKAWAVEVYKQARPGYHGLAQGTVDGVLK'


    macromolecule = MacroMol(['> proteina', prot], orientAmino, [0,2,3,4,6,11,23], 0, False, True)
    print(macromolecule.result)

    sequences = [['> a ',
                  'ccgcgcgtactcagtactgcagtcagtaCGGTGTGTGTGTGTGTgacgctacgCCCCCAAAAAatcagcatcagactcagcatagctcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctatatatatgctacgtacgcatgactcagactaccg'], \
                 ['>b ',
                  'ccgcgcgtactcCGTGTGTagtactgcagtcagtagtcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctatatatatgctacgtacgcatgactcagactaccg'],
                 ['> c ',
                  'ccgcgcgtactcagtactgcagtcagtagtcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctacgatcgatcgatcgatcgatcgatcgatcatatatatgctacgtacgcatgactcagactaccg']]

    c=process_all(sequences,orientaseq,[0,1,2,2],2)
    print(c.__len__())
    print(c)
    print('-----------------------------------------')
    sequences = readfilefasta('ITS1segment_C.fasta')
    print(process_all(sequences, orientaseq, [0], 2))





