#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datapond
import math
import numpy as np
from multiprocessing import Pool


class allfromctrecord():
    def __init__(self, filenameseq=''):
        self.filename = filenameseq
        self.ct = []
        self.connectivity = []
        self.sequence = ''
        self.name = ''
        self.n = 0
        self.openctfile()

    def add(self):
        self.ct.append([self.name, self.sequence, self.connectivity])

    def openctfile(self):
        f = open(self.filename, 'r')
        # print ('==>', self.filename)
        # splitting firstline to get sequence size and name
        counter = 0
        for line in f:
            if 'dG' in line:
                counter += 1
                self.connectivity = []
                firstline = line.split('dG')
                self.n = int(firstline[0])
                # print self.n
                self.name = firstline[1].strip()[12:] + '_' + str(counter)
                self.sequence = ''
                cnumber = 0

            else:
                currentline = line
                cnumber = int(currentline[:7])
                nucleotide = currentline[7]
                self.sequence = self.sequence + nucleotide
                previousN = int(currentline[8:15])
                nextN = int(currentline[15:23])
                bondedto = int(currentline[23:30])
                # print previousN,nextN,bondedto
                if nextN != 0:
                    self.connectivity.append((previousN, nextN - 1))
                if (bondedto != 0) and (bondedto - 1 > previousN):
                    self.connectivity.append((previousN, bondedto - 1))

            if cnumber == self.n:
                self.add()


class ctrecord():
    def __init__(self, filenameseq=''):
        self.filename = filenameseq
        self.connectivity = []
        self.results = []
        self.sequence = ''
        self.ponderations = {}
        self.name = 'test'  # sequence name
        self.n = 0  # sequence size
        if filenameseq != '':
            self.openctfile()  # file to process

    def myInitCalc(self, ponderations, nameseq, sequence, connectivity):
        self.outdata = []
        self.name = nameseq
        self.sequence = sequence.upper()
        self.connectivity = connectivity
        self.n = len(self.sequence)
        self.ponderations = ponderations
        self.calc()
        s = self.out_header()
        self.outdata.append(s)
        self.outdata.append(self.results)
        return self.outdata

    def write_fasta_format(self):
        return '>' + self.name + '\n' + self.sequence.upper() + '\n'

    def __str__(self):
        '''
        return sequence into fasta format
        '''
        return '>' + self.name + '\n' + self.sequence.upper() + '\n'

    def setponderations(self, ponderations):
        '''
        set ponderation dictionary
        '''
        self.ponderations = ponderations

    def getseqSize(self):
        '''
        return sequence size
        '''
        return self.n

    def getconnections(self):
        '''
        return connectivity
        '''
        return self.connectivity

    def openctfile(self):
        '''
        get sequence and connectivity from ct file
        '''
        f = open(self.filename, 'r')
        # splitting firstline to get sequence size and name
        firstline = f.readline().split('dG')
        self.n = int(firstline[0])
        self.name = firstline[1].strip()[12:]
        cnumber = 0
        while cnumber != self.n:
            currentline = f.readline()
            cnumber = int(currentline[:7])
            nucleotide = currentline[7]
            self.sequence = self.sequence + nucleotide
            previousN = int(currentline[8:15])
            nextN = int(currentline[15:23])
            bondedto = int(currentline[23:30])
            # print previousN,nextN,bondedto
            if nextN != 0:
                self.connectivity.append((previousN, nextN - 1))
            if (bondedto != 0) and (bondedto - 1 > previousN):
                self.connectivity.append((previousN, bondedto - 1))

                # np.savetxt('salidaconectividadx.csv', self.connectivity, delimiter=',')

    def setnameSeq(self, nameseq):
        self.name = nameseq

    def builtmatrix(self):

        self.Matrix = np.zeros([len(self.connectivity), len(self.connectivity)])
        # print '------->',len(self.connectivity), type(self.connectivity)
        # print list(''.join(self.connectivity))
        for i in range(len(self.connectivity)):
            for j in range(len(self.connectivity)):
                if j != i:
                    x = self.connectivity[i]
                    y = self.connectivity[j]

                    if x[0] == y[0] or x[0] == y[1] or x[1] == y[0] or x[1] == y[1]:
                        self.Matrix[i, j] = self.Matrix[j, i] = 1.0

        for i in range(len(self.connectivity)):
            self.Matrix[i, i] = 0.0

        '''
        np.savetxt('salida.csv', self.Matrix,delimiter=',')
        '''

        counter = {}
        for i in range(len(self.sequence)):
            counter[i] = 0

        matrix2 = np.zeros([len(self.sequence), len(self.sequence)])

        for x in self.connectivity:
            counter[x[0]] += 1
            counter[x[1]] += 1
            matrix2[x[0], x[1]] = 1.0

        # c_INDEX
        c = list(counter.values())

        # print 'valores a evaluar', len(c)

        sum = 0.0
        # print '------->>>>>', len(matrix2)

        for i in range(len(matrix2)):
            for j in range(len(matrix2)):
                if (j > i) and (matrix2[i, j] != 0.0):
                    sum = sum + 1.0 / math.sqrt(c[i] * c[j])

        self.results.append(self.name)
        self.results.append(len(self.connectivity))
        self.results.append(sum)
        self.results.append(self.ponderations['name'])

    def toponderate_doit(self):
        Edgespond = []
        for i in range(len(self.connectivity)):
            x = self.connectivity[i][0]
            y = self.connectivity[i][1]

            pond = self.ponderations[self.sequence[x]] + self.ponderations[self.sequence[y]]
            # print(self.sequence[x], self.sequence[y],pond)
            Edgespond.append(pond)
        # print(self.ponderations)
        for i in range(len(self.Matrix)):
            for j in range(len(self.Matrix)):
                if j > i:
                    if self.Matrix[i, j] == 1.0:
                        self.Matrix[i, j] = self.Matrix[j, i] = (Edgespond[i] + Edgespond[j]) / 2.0
                if i == j:
                    self.Matrix[i, i] = Edgespond[i]

        sumcol = np.zeros(len(self.Matrix))
        for i in range(len(self.Matrix)):
            sumcol[i] = self.Matrix[i][:].sum()
            if sumcol[i] == 0:
                sumcol[i] = 1.0

        self.Matrix = self.Matrix / sumcol
        # print self.Matrix
        auxmatrix = self.Matrix.copy()
        self.results.append(round(self.Matrix.trace(), 4))

        for i in range(14):
            auxmatrix = auxmatrix.dot(self.Matrix)
            self.results.append(round(auxmatrix.trace(), 4))

    def calc(self):
        self.builtmatrix()
        self.toponderate_doit()

    def out_header(self):
        s = "Name" + "\t" + "Edge Number" + "\t" + "Edge connnectivity Index" + "\t" + "Ponderation" + "\t"
        for i in range(1, 16):
            s = s + "Momment" + str(i) + "\t"
        s = s + '\n'
        return s

    def saveout(self, fileoutput):
        f = open(fileoutput, 'w')
        s = self.out_header()
        f.write(s)
        f.write(str(self.results))

        f.close()


class ct2record(ctrecord):
    def __init__(self, ponderations, nameseq, sequence, connectivity):
        ctrecord.__init__(self, filenameseq='')
        self.outdata = []
        self.name = nameseq
        self.sequence = sequence.upper()
        self.connectivity = connectivity
        self.n = len(self.sequence)
        self.ponderations = ponderations
        self.calc()
        s = self.out_header()
        self.outdata.append(s)
        self.outdata.append(self.results)

    def getresults(self):
        return self.outdata


def calculator(sequence, ponderation):
    if isinstance(sequence['connection'], str):
        s = sequence['connection'].replace('), (', ';')
        s = s.replace('[(', ' ')
        s = s.replace(')]', ' ')
        s = s.replace('[', '')
        s = s.replace(']', '')
        allconect = s.split(';')
        connections = []
        for x in allconect:
            data = (x.split(','))

            x1 = int(data[0])
            x2 = int(data[1])
            connections.append((x1, x2))
    else:
        connections = sequence['connection']

    dna2mol = ct2record(ponderation, sequence['name'], str(sequence['seq']), connections)

    return dna2mol.results


def process_all(sequences, ponderations):
    print(ponderations, 'I was here already..how change ponderation.')
    temp = []
    for pk in ponderations:
        pondname, data = datapond.getPond(pk, 2)
        print(pondname, data)
        ponderation = {'name': pondname, 'A': data[0], 'C': data[3], 'T': data[2], 'G': data[1], 'U': data[4],
                       'W': 0.20, 'N': 0.20}
        pool = Pool()
        resultsx = [pool.apply_async(calculator, (x, ponderation,), ) for x in sequences]
        [temp.append(r.get()) for r in resultsx]

    return temp


if __name__ == '__main__':
    '''
    rnamol = ctrecord('/home/reymolina/PycharmProjects/calc/3AAC000003.ct')
    #rnamol = ctrecord('/home/reymolina/PycharmProjects/calc/Ctfiles/ITS2segment_TKT_pure_dep.fasta.ct')
    #print rnamol
    myown = {'name':'my_charges', 'A':0.22, 'C':0.19,'T':0.21, 'G':0.24, 'U':0.21 }
    #print rnamol.getseqSize()
    #print rnamol.getconnections()
    rnamol.setponderations(myown)
    rnamol.calc()
    print rnamol.results
    rnamol.out_header()
    rnamol.saveout('operationcomplete.csv')
    '''

    # direct calculation
    myown = {'name': 'my_charges', 'A': 0.22, 'C': 0.19, 'T': 0.21, 'G': 0.24, 'U': 0.21}
    sequence = 'cgacgtcatgcagtacgtca'
    print(len(sequence))
    connection = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 15), (7, 8), (8, 9), (8, 16), (9, 10),
                  (10, 11), (11, 12),
                  (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19)]
    dna2mol = ct2record(myown, 'test', sequence, connection)
    print(dna2mol.getresults())
