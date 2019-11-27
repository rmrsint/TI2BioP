#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'reymolina'


import re
import sys
import numpy as np
from multiprocessing import cpu_count, Process, Pool, Value, Array

def partez(bct):

    istart = []  # stack of indices of opening parentheses
    d = {}
    bonds = bct

    for i, c in enumerate(bonds):
        if c == '(':
            istart.append(i)
        if c == ')':
            try:
                d[istart.pop()] = i
            except IndexError:
                print('Too many closing parentheses')
    if istart:  # check if stack is empty afterwards
        print('Too many opening parentheses')

    bonded = [(x, d[x]) for x in d.keys()]

    return bonded

def parenthesesmatch(bct, sizeseq):
    istart = []  # stack of indices of opening parentheses
    d = {}
    bonds = bct
    bondx = []
    for i in range(sizeseq-1):
        bondx.append((i,i+1))


    for i, c in enumerate(bonds):
        if c == '(':
            istart.append(i)
        if c == ')':
            try:
                d[istart.pop()] = i
            except IndexError:
                print('Too many closing parentheses')
    if istart:  # check if stack is empty afterwards
        print('Too many opening parentheses')

    bonded = bondx + [(x, d[x]) for x in d.keys()]
    print (bonded)
    return bonded



class allfromxfasta():
    def __init__(self, filenameseq):
        self.filename = filenameseq
        self.xfasta = []
        self.openxfasta()
        self.processall()


    def processall(self):
        conectz = []
        for i in range((len(self.xfasta))):
            conectz.append(self.xfasta[i][2])
        res = Pool().map(partez, conectz)

        for i in range(len(res)):
            bonded = []
            for j in range (len(self.xfasta[i][2])-1):
                bonded.append((j, j+1))

            self.xfasta[i][2] = bonded + res[i]


    def openxfasta(self):
        f = open(self.filename, 'r')
        for line in f:
            if re.search('(\>.+)', line):
                label = line.split(' ')[0][1:]

            elif re.search('(^[A*a*T*t*C*c*G*g*U*u*]*$)', line):
                seq = line.rstrip()

            elif '(' in line:
                splitline = line.split(' ')
                energy = float(re.findall(r'\(([\d\.-]+)\)', splitline[1])[0])
                conectivity = splitline[0].rstrip()
                self.xfasta.append([label, seq, conectivity , energy, ])

        f.close()

