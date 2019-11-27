__author__ = 'reymolina'
'''
Calculation of topological indexes for BioPolymers
Four Colors Map
'''

import Bio


from multiprocessing import Pool

from pylab import *
# from Drawspiral import MyFCMFrame
from Drawsp_qt import MyWidget

import datapond


Aminoup = 'AVLIPFWM'  # red
AminoDown = 'GSTCYNQ'  # Green
AminoLeft = 'DEBZX'  # Yellow
AminoRight = 'KRH'  # Blue


def checktouch(cluster00, cluster01):
    contact = 0
    cluster1 = cluster00[0]
    cluster2 = list(cluster01[0])

    for elem in cluster1:
        row, col = elem
        if (row + 1, col) in cluster2:
            contact = 1
        if (row - 1, col) in cluster2:
            contact = 1
        if (row, col + 1) in cluster2:
            contact = 1
        if (row, col - 1) in cluster2:
            contact = 1
    return contact


def matrixnumber(n):
    if n != 0:
        root = int(math.sqrt(n))
        if root ** 2 != n:
            root += 1
        current = root ** 2
        # print ('secuence length  = ', n)
        # print ('current ', current)
        # print ('root ', root)
        cell_list = np.array([])
        if root % 2 == 0:
            colx = root - 1
            cell_list = np.append(cell_list, {'row': root - 1, 'col': root - 1, 'pos': current})

            for rowx in range(root - 2, -1, -1):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
            rowx = 0
            for colx in range(root - 2, -1, -1):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
            colx = 0
            for rowx in range(1, root):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
            rowx = root - 1
            for colx in range(1, root - 1):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
        else:
            cell_list = np.append(cell_list, {'row': 0, 'col': 0, 'pos': current})
            colx = 0
            for rowx in range(1, root):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
            rowx = root - 1
            for colx in range(1, root):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
            colx = root - 1
            for rowx in range(root - 2, -1, -1):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
            rowx = 0
            for colx in range(root - 2, 0, -1):
                current = current - 1
                cell_list = np.append(cell_list, {'row': rowx, 'col': colx, 'pos': current})
    return cell_list, root


def totaladj_index(mcontactos):
    counter = 0
    for i in range(len(mcontactos)):
        for j in range(len(mcontactos)):
            if j > i and mcontactos[i][j] == 1:
                counter += 1
    return 2 * counter


class MacroMol():
    def __init__(self, seqraw, ponderations,  seqsort=0, calc=True):
        self.sequence = ''
        self.n = 0
        self.names = seqraw[0]
        self.seqarray = []
        self.topindexs = []
        self.macrotyp = seqsort
        self.ponderation = []
        self.ponds = ponderations
        self.spiral = []
        self.fourcolormap = []
        self.cluster_yellow_list = []
        self.cluster_red_list = []
        self.cluster_green_list = []
        self.cluster_blue_list = []
        self.result = []
        self.result.append(self.names)
        self.all_in_one(seqraw[1],  seqsort, perfcalc=calc)


    def readfromfile(self, filename):
        f1 = open(filename, 'r')
        self.seqclear = []
        self.seqarray = f1.readlines()
        f1.close()

        for x in self.seqarray:
            if '>' not in x:
                self.seqclear.append(x)


    def getseqat(self, n):
        self.names = self.seqarray[2 * n]
        self.sequence = self.seqarray[2 * n + 2]
        self.n = len(self.sequence)

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
        self.setClusters()
        self.n = self.sequence.__len__()

    def listcontact(self, nx):
        counter = -1
        for i in range(len(self.connectivity)):
            counter += self.connectivity[nx][i]
        return counter

    def setClusters(self):
        if (self.macrotyp == 0 ):
            self.cluster_yellow = AminoLeft
            self.cluster_green = AminoDown
            self.cluster_blue = AminoRight
            self.cluster_red = Aminoup
        elif self.macrotyp == 1:
            self.cluster_yellow = 'A'
            self.cluster_green = 'C'
            self.cluster_blue = 'T'
            self.cluster_red = 'G'
        elif self.macrotyp == 2:
            self.cluster_yellow = 'A'
            self.cluster_green = 'C'
            self.cluster_blue = 'U'
            self.cluster_red = 'G'

        '''
          macrotyp: 0 proteins
          macrotyp: 1 DNA
          macrotyp: 2 RNA

        '''
    def setPond(self, data, printout=False):

        clusternumber = self.matrixcluster.max() + 1

        if printout:
            np.savetxt('vermatrix.txt', self.matrixcluster, fmt='%2.0f')
            f = open('verquepasa.txt', 'w')
            # print self.gridlabel

            for i in self.gridlabel:
                f.write(str(i))
                

            f.close()
            print(len(self.gridlabel))
            print(len(self.gridlabel[0]))
            print(len(self.matrixcluster))

        ponderation =[]
        for i in range(clusternumber):
            ponderation.append(0)

        for i in range(len(self.matrixcluster)):
            for j in range(len(self.matrixcluster)):
                if self.matrixcluster[i][j] != -1:

                    try:
                        if self.macrotyp == 0:
                            ponderation[self.matrixcluster[i][j]] += data[datapond.AASEQ2.index(self.gridlabel[i][j][0])]
                        else:
                            ponderation[self.matrixcluster[i][j]] += data[datapond.NAcids.index(self.gridlabel[i][j][0])]
                    except:
                        print(i, j, self.gridlabel[i][j][0])
                        print(len(self.matrixcluster))


        return ponderation


    def setcalc(self, n):
        self.macrotyp = n
        self.macrotyp = n

    def showspiral(self):
        figure(1)
        imshow(self.fourcolormap, interpolation='nearest')  # nearest
        # imshow(self.fourcolormap, aspect='auto', cmap=plt.get_cmap('gist_ncar'))
        grid(True)
        show()

    def showspiral2(self):
        presentar = MyWidget(title='Four colors map representation'+'  '+self.names, savefilename =self.names + '_fcm_graph.png', spiralx=self.fourcolormap, labels=self.gridlabel)
        print(self.names)
        presentar.show()



    def showmol(self):
        print(self.sequence)
        print(self.macrotyp)
        print('Sequence length ==> ', self.n)
        print('Blue cluster ==>', self.cluster_blue)
        print('Red cluster ==>', self.cluster_red)
        print('Yellow cluster ==>', self.cluster_yellow)
        print('Green cluster ==>', self.cluster_green)

    def showcluster(self):
        print('Blue cluster -->', len(self.cluster_blue_list))
        print(self.cluster_blue_list)
        print('Yellow Cluster -->', len(self.cluster_yellow_list))
        print(self.cluster_yellow_list)
        print('Red Cluster -->', len(self.cluster_red_list))
        print(self.cluster_red_list)
        print('Green Cluster -->', len(self.cluster_green_list))
        print(self.cluster_green_list)

    def clustering(self):

        for kcluster in range(4):

            if kcluster == 0:
                ktproc = self.cluster_blue_list
                size = self.cluster_blue_list.__len__()
            elif kcluster == 1:
                ktproc = self.cluster_yellow_list
                size = self.cluster_yellow_list.__len__()
            elif kcluster == 2:
                ktproc = self.cluster_green_list
                size = self.cluster_green_list.__len__()
            elif kcluster == 3:
                ktproc = self.cluster_red_list
                size = self.cluster_red_list.__len__()

            # build a list o empty clusters
            newdistr = []
            if size != 0:

                for i in range(size):
                    newdistr.append([])

                lastelement = False
                I = 0

                while not lastelement:

                    row, col = ktproc[0]
                    ktproc.remove(ktproc[0])
                    size = len(ktproc)
                    newdistr[I].append((row, col))

                    if (row, col - 1) in ktproc:
                        newdistr[I].append((row, col - 1))
                        ktproc.remove((row, col - 1))
                        size = len(ktproc)

                    if (row - 1, col) in ktproc:
                        newdistr[I].append((row - 1, col))
                        ktproc.remove((row - 1, col))
                        size = len(ktproc)

                    if (row, col + 1) in ktproc:
                        newdistr[I].append((row, col + 1))
                        ktproc.remove((row, col + 1))
                        size = len(ktproc)

                    if (row + 1, col) in ktproc:
                        newdistr[I].append((row + 1, col))
                        ktproc.remove((row + 1, col))
                        size = len(ktproc)

                    I += 1

                    if len(ktproc) == 0:
                        lastelement = True

                # print 'Before -->', len(newdistr)
                while [] in newdistr:
                    newdistr.remove([])

                # Rebuilt clusters

                for i in range(len(newdistr)):

                    for x in newdistr[i]:

                        y, z = x

                        for j in range(i + 1, len(newdistr)):
                            if (y, z + 1) in newdistr[j]:
                                newdistr[i].append((y, z + 1))
                                newdistr[j].remove((y, z + 1))

                            if (y, z - 1) in newdistr[j]:
                                newdistr[i].append((y, z - 1))
                                newdistr[j].remove((y, z - 1))

                            if (y + 1, z) in newdistr[j]:
                                newdistr[i].append((y + 1, z))
                                newdistr[j].remove((y + 1, z))

                            if (y - 1, z) in newdistr[j]:
                                newdistr[i].append((y - 1, z))
                                newdistr[j].remove((y - 1, z))

                while [] in newdistr:
                    newdistr.remove([])

            if kcluster == 0:
                self.cluster_blue_list = newdistr
            elif kcluster == 1:
                self.cluster_yellow_list = newdistr
            elif kcluster == 2:
                self.cluster_green_list = newdistr
            elif kcluster == 3:
                self.cluster_red_list = newdistr

    def updatematrixcluster(self):
        matrixsize = self.cluster_red_list.__len__() + self.cluster_green_list.__len__() + self.cluster_blue_list.__len__() + self.cluster_yellow_list.__len__()
        if self.macrotyp == 0:
            self.result.append(matrixsize)
        else:
            self.result.append(matrixsize * 3)

        self.connectivity = np.zeros([matrixsize, matrixsize], dtype=int)
        bigmatrix = []
        for i in range(matrixsize):
            bigmatrix.append([])


        # print 'number of cluster ', matrixsize
        root = int(math.sqrt(len(self.sequence)))
        if root ** 2 != len(self.sequence):
            root += 1

        self.matrixcluster = np.zeros([root, root], dtype=int)
        self.matrixcluster -= 1

        N = 0

        for i in self.cluster_red_list:
            bigmatrix[N].append(i)
            for elem in i:
                row, col = elem
                self.matrixcluster[row][col] = N
            N += 1

        for i in self.cluster_green_list:
            bigmatrix[N].append(i)
            for elem in i:
                row, col = elem
                self.matrixcluster[row][col] = N
            N += 1

        for i in self.cluster_blue_list:
            bigmatrix[N].append(i)
            for elem in i:
                row, col = elem
                self.matrixcluster[row][col] = N
            N += 1

        for i in self.cluster_yellow_list:
            bigmatrix[N].append(i)
            for elem in i:
                row, col = elem
                self.matrixcluster[row][col] = N
            N += 1

        # print self.matrixcluster

        for i in range(matrixsize):
            for j in range(matrixsize):
                if i == j:
                    self.connectivity[i][j] = 1
                elif j > i:
                    self.connectivity[i][j] = checktouch(bigmatrix[i], bigmatrix[j])
                    self.connectivity[j][i] = self.connectivity[i][j]

                    # print(self.connectivity)
                    #np.savetxt('salida.txt', self.connectivity, fmt = '%2.0f')


    def definepond(self):
        # aminoacids one letter code
        if self.macrotyp == 0:
            self.idstr = Bio.Data.IUPACData.protein_letters
        elif self.macrotyp == 1:
            # DNA
            self.idstr = 'CGAT'
        elif self.macrotyp == 2:
            # RNA
            self.idstr = 'CGAU'

    def buildpondcluster(self, printout=False):
        clusternumber = self.matrixcluster.max() + 1

        if printout:
            np.savetxt('vermatrix.txt', self.matrixcluster, fmt='%2.0f')
            f = open('verquepasa.txt', 'w')
            # print self.gridlabel

            for i in self.gridlabel:
                f.write(str(i))

            f.close()
            print( len(self.gridlabel))
            print( len(self.gridlabel[0]))
            print( len(self.matrixcluster))

        for i in range(clusternumber):
            self.ponderation.append(0)

        for i in range(len(self.matrixcluster)):
            for j in range(len(self.matrixcluster)):
                if self.matrixcluster[i][j] != -1:
                    # print self.gridlabel[i][j][0]
                    try:
                        if self.macrotyp == 0:
                            self.ponderation[self.matrixcluster[i][j]] += datapond.amino_acids[self.gridlabel[i][j][0]]
                        else:
                            self.ponderation[self.matrixcluster[i][j]] += datapond.nitrog_bases[self.gridlabel[i][j][0]]
                    except:
                        print(i, j)
                        print(len(self.gridlabel))
                        print(len(self.matrixcluster))

                        # print self.ponderationn

    def buildspiral(self):
        StepW = 0
        NN = self.n
        # print NN
        cell_list, root = matrixnumber(NN)
        self.grid = np.zeros([root, root])
        self.gridlabel = []
        while (NN > 0) and (NN != 0):
            for i in cell_list:
                self.grid[int(i['row']) + StepW][StepW + int(i['col'])] = int(i['pos'])
                NN = int(i['pos']) - 1
            if NN != 0:
                cell_list, root = matrixnumber(NN)
            StepW += 1
        root = int(math.sqrt(self.n))
        if root ** 2 < self.n:
            root += 1
        self.fourcolormap = np.zeros([root, root], dtype=int)

        for i in range(len(self.grid)):
            c = []
            for j in range(len(self.grid)):
                if self.grid[i][j] - 1 < self.n:
                    elem = int(self.grid[i][j] - 1)
                    simbol = self.sequence[elem]
                    c.append(simbol + str(elem + 1))
                    self.fourcolormap[i][j] = -1
                    if simbol in self.cluster_blue:
                        self.fourcolormap[i][j] = 1
                        self.cluster_blue_list.append((i, j))
                    elif simbol in self.cluster_green:
                        self.fourcolormap[i][j] = 4
                        self.cluster_green_list.append((i, j))
                    elif simbol in self.cluster_red:
                        self.fourcolormap[i][j] = 2  # 0.9227
                        self.cluster_red_list.append((i, j))
                    elif simbol in self.cluster_yellow:
                        self.fourcolormap[i][j] = 3
                        self.cluster_yellow_list.append((i, j))
            self.gridlabel.append(c)

        # fixing labels
        # print ' root = ', root

        c1 = self.gridlabel[0]
        c2 = self.gridlabel[-1]

        if int(c1[-1][1:]) > int(c1[0][1:]):
            maxc1 = int(c1[-1][1:])
        else:
            maxc1 = int(c1[0][1:])

        if int(c2[-1][1:]) > int(c2[0][1:]):
            maxc2 = int(c2[-1][1:])
        else:
            maxc2 = int(c2[0][1:])

        if len(c1) < root and len(c2) < root:

            for i in range(len(self.gridlabel)):
                if len(self.gridlabel[i]) < root:
                    for j in range(root - len(self.gridlabel[i])):
                        if maxc1 > maxc2:
                            self.gridlabel[i].append('--')
                        else:
                            self.gridlabel[i].insert(0, '--')

        if len(c1) == root and len(c2) == root:
            for i in range(len(self.gridlabel)):
                if len(self.gridlabel[i]) < root:
                    for j in range(root - len(self.gridlabel[i])):
                        self.gridlabel[i].append('--')

        if len(c1) == root and len(c2) < root:
            for i in range(len(self.gridlabel)):
                if len(self.gridlabel[i]) < root:
                    for j in range(root - len(self.gridlabel[i])):
                        self.gridlabel[i].append('--')

        if len(c1) < root and len(c2) == root:

            for i in range(len(self.gridlabel)):
                if len(self.gridlabel[i]) < root:
                    for j in range(root - len(self.gridlabel[i])):
                        self.gridlabel[i].insert(0, '--')

                        # pprint.pprint(self.gridlabel)
                        # print self.gridlabel

    def calc(self):
        totaladj = totaladj_index(self.connectivity)

        for k in self.ponds:
            #print k, 'sadsd'
            result = []
            pondname,data = datapond.getPond(k, self.macrotyp)

            result.append(totaladj)
            result.append(pondname)
            self.ponderation = self.setPond(data)


            matrixcalc = np.zeros([len(self.connectivity), len(self.connectivity)])

            for i in range(len(self.connectivity)):
                for j in range(len(self.connectivity)):
                    if j >= i:
                        if self.connectivity[i][j] == 1:
                            matrixcalc[i][j] = (self.ponderation[i] + self.ponderation[j] ) / 2
                            matrixcalc[j][i] = matrixcalc[i][j]
            sumcol = np.zeros(len(matrixcalc))

            for i in range(len(matrixcalc)):
                sumcol[i] = matrixcalc[i][:].sum()
                if sumcol[i] == 0.0:
                    sumcol[i] = 1.0

            matrixcalc = matrixcalc / sumcol
            auxmatrix = matrixcalc.copy()
            result.append(round(matrixcalc.trace(), 4))

            for i in range(14):
                auxmatrix = auxmatrix.dot(matrixcalc)
                result.append(round(auxmatrix.trace(), 4))

            self.result.append(result)

    def all_in_one(self, sequence, ksort, perfcalc=False):
        self.addsequence(sequence, ksort)
        # print self.n
        self.buildspiral()
        self.clustering()
        self.updatematrixcluster()
        self.buildpondcluster()
        if perfcalc:
            self.calc()

        # print 'RESULTS.....>', self.result


def save_results(prefixname, results):
    print(results.__len__())



def calcdescriptors(seqdata,ponderations,  seqsort):
    mol = MacroMol(seqdata, ponderations, seqsort,  calc=True)
    return mol.result

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



def process_all(sequences, ponderations, seqsort):
    pool = Pool()
    temp = []
    """
    print('-------------------------------')
    print(ponderations)
    print('-------------------------------')
    print(seqsort)
    print('-------------------------------')
    print(sequences)
    print('-------------------------------')
    print('I was )    
    """

    resultx = [pool.apply_async(calcdescriptors, args=(x, ponderations, seqsort, )) for x in sequences]
    [temp.append(r.get()) for r in resultx]

    return temp

if __name__ == "__main__":
    macromolecule = MacroMol(['> xx',
                              'ccgcgcgtactcagtactgcagtcagtagtcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctatatatatgctacgtacgcatgactcagactaccg'],
                             [0],2, True)
    print(macromolecule.result)
    macromolecule.showspiral2()
    sequences = [['> a ',
                  'ccgcgcgtactcagtactgcagtcagtaCGGTGTGTGTGTGTGTgacgctacgCCCCCAAAAAatcagcatcagactcagcatagctcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctatatatatgctacgtacgcatgactcagactaccg'], \
                 ['>b ',
                  'ccgcgcgtactcCGTGTGTagtactgcagtcagtagtcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctatatatatgctacgtacgcatgactcagactaccg'],
                 ['> c ',
                  'ccgcgcgtactcagtactgcagtcagtagtcaggccgctctgctgccttttttttctttggcgtgcgtgcgtgcgacacgtcgctcgatcgatcgactgatgctcgtacgtagatcgactgctacgatcgatcgatcgatcgatcgatcgatcatatatatgctacgtacgcatgactcagactaccg']]

    print(process_all(sequences,[0],2))

    sequences = readfilefasta('ITS1segment_C.fasta')
    results=process_all(sequences,[0],2)
    save_results('salida', results)
    #print results

    #sequences = readfilefasta('ITS1segment_C.fasta')
    #print process_all(sequences,2)

    #sequences = readfilefasta('ITS1segment_C.fasta')
    #print process_all(sequences,2)



    # macromolecule.showspiral2()
    #macromolecule.showspiral2()








