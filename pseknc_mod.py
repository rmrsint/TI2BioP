#!/usr/bin/python

# -----------------------------------------
# PseKNC - 06/11/2014
# Jordan Brooker
# Jrdnbrkr@gmail.com
# -----------------------------------------

from __future__ import division
import Data
import math, string, re, sys, getopt

#  _mean: Calculates mean value of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the mean value of listy

def _mean(listy):
    return float(sum(listy)) / len(listy)


#  _std(listy): Calculates standard deviation of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the standard deviation of listy

def _std(listy, ddof=1):
    mean = _mean(listy)
    # for i in listy, (i - mean)^2 => list
    temp = [math.pow(i - mean, 2) for i in listy]
    # temp is added together, divided by length of listy,
    # and the square root is taken
    res = math.sqrt(float(sum(temp)) / len(listy))
    return res


#  sepSequence
# -----------------------------------------
# inputs: seq, string and k, int
# output: list of k-tuples

def sepSequence(seq, k):
    i = k - 1
    seqq = []
    while i < len(seq):
        j = 0
        nuc = ''
        while j < k:
            nuc = seq[i - j] + nuc
            j = j + 1
        seqq.append(nuc)
        i += 1
    return seqq


#  getValues: Returns a line of values for one property
#             from physicochemical property files
# -----------------------------------------
# input: prop = string of one property and supInfo = string
# output: values = a string representing a list of property values

def getValues(prop, supInfo):
    values = ""
    name = re.search(prop, supInfo)
    if name:
        strr = prop + '\s*\,(.+)'
        b = re.search(strr, supInfo)
        if b:
            values = b.group(1)
    return values


#  getSpecificValue: Returns a property value for a specific di- or tri-
#                    nucleotide
# -----------------------------------------
# input: olinuc = string, prop = string, supInfo = string
# output: value = an int that is property value of an olinuc

def getSpecificValue(olinuc, olinucs, prop, values, supInfo):
    values = values.split(",")
    #valueS = [float(x) for x in values.split(",")]
    count = olinucs.index(olinuc)
    value = values[count]
    return float(value)


#  hn: Hn function
# -----------------------------------------
# inputs: olinuc = string, prop = string, supInfo = string
# output: temp = int

def hn(olinuc, olinucs, prop, supInfo, values):
    #values = getValues(prop,supInfo).rstrip()
    h0 = float(getSpecificValue(olinuc, olinucs, prop, values, supInfo))
    valueS = [float(x) for x in values.split(",")]
    temp = float((h0 - _mean(valueS)) / _std(valueS))
    return temp


#  theta2: Theta(i,i+j) function, for Type I
# -----------------------------------------
# input: seq = string, props = string, i = int,
#	 Type = int, k = int, j = int, supInfo = string
# output: summ = int

def theta2(seq, olinucs, props, i, k, j, supInfo):
    summ = 0
    values = ''
    for prop in props:
        values = getValues(prop, supInfo).rstrip()
        hn1 = hn(seq[i - 1], olinucs, prop, supInfo, values)
        hn2 = hn(seq[i + j - 1], olinucs, prop, supInfo, values)
        subsqr = math.pow(float(hn1 - hn2), 2)
        summ = summ + subsqr
    return float(summ) / len(props)


#  J: J(i,i+j) function, for Type II
# --------------------------------------
# inputs: seq = string, prop = string, i = int, 
#         k = int, j = int, supInfo = string
# output: product = int

def J(seq, olinucs, prop, i, k, j, supInfo):
    values = getValues(prop, supInfo)
    hn1 = hn(seq[i - 1], olinucs, prop, supInfo, values)
    hn2 = hn(seq[i + j - 1], olinucs, prop, supInfo, values)
    return float(hn1 * hn2)


#  theta1: Theta(j) and Tau(LamGam) function
# -----------------------------------------
# input: seq = string, props = string, Type = int
#        k = int, j = int, supInfo = string
# output: final = int

def theta1(seq, olinucs, props, Type, k, j, supInfo):
    k = int(k)
    gamma = len(props)
    seqq = sepSequence(seq, k)
    i = 1
    a = 0
    if Type == 1:
        var = len(seq) - int(k) - int(j) + 1
        while i <= var:
            b = 0
            b = theta2(seqq, olinucs, props, i, k, j, supInfo)
            a = a + b
            i = i + 1
    else:
        ii = 0
        var = len(seq) - int(k) - int(j / gamma)
        while i <= var:
            b = 0
            if ii == gamma:
                ii = 0
                b = J(seqq, olinucs, props[ii], i, k, int(j / gamma), supInfo)
                a = a + b
            else:
                b = J(seqq, olinucs, props[ii], i, k, int(j / gamma), supInfo)
                a = a + b
            ii = ii + 1
            i = i + 1
    final = float(a) / var
    return final


#  pseKNCHelper: Creates list of adjusted frequency values for 4^k
#  oligonucleotides and the lambda terms
# -----------------------------------------
# input: seq = string, props = string, Type = int, k = int, j = int, w = int
# output: freqs = list of ints

def pseKNCHelper(seq, olinucs, props, Type, k, j, w, geneticMaterial, supInfo):
    gamma = len(props)
    freqs = []
    seqq = sepSequence(seq, k)
    olinucs = olinucs.split(",")
    for olinuc in olinucs:
        freq = seqq.count(olinuc)
        freq = float(freq / len(seqq))
        total = 0
        i = 1
        if Type == 2:
            j = j / gamma
        while i <= j:
            total = total + theta1(seq, olinucs, props, Type, k, i, supInfo)
            i = i + 1
        total = float(freq / (1 + (float(w) * total)))
        total = int(total * 1000) / 1000.0
        freqs.append(total)
    #Computing Lambda terms...
    fourK = math.pow(4, k)
    mu = fourK + 1
    while (fourK + 1) <= mu <= (fourK + j):
        top = float(w) * theta1(seq, olinucs, props, Type, k, int(mu - fourK), supInfo)
        bottomTheta = 0
        bottom = 0
        i = 1
        while 1 <= i <= j:
            bottomTheta += theta1(seq, olinucs, props, Type, k, i, supInfo)
            i += 1
        bottom = 1 + (float(w) * bottomTheta)
        term = float(top / bottom)
        term = int(term * 1000) / 1000.0
        freqs.append(term)
        mu += 1
    return freqs


#  pseKNC: Opens input and output files, calls functions to calculate
#  values and writes outputs to the output file in appropriate format
# -----------------------------------------
# input: inputFile = string, outputFile = string, propNames = string,
#        Type = int, k = int, j = int, w = int, formatt = string
# output: nothing, write to output file

def pseKNC(sequenceX, props, Type, k, j, w, formatt, geneticMaterial):
    j = float(j)
    listy = ''
    output = ''
    outputs = ''
    gamma = len(props)
    #Getting supporting info from files
    if k == 2:
        if geneticMaterial == 'DNA':
            supInfo = Data.SetS1DNAInfo
        else:
            supInfo = Data.SetS1RNAInfo
    else:
        supInfo = Data.SetS3DNAInfo

    #SupFile = open(SupFileName,'r')
    #supInfo = SupFile.read()
    o = re.search('Physicochemical properties\,(.+)\n', supInfo)

    olinucs = ''
    if o:
        olinucs = o.group(1).rstrip()

    if Type == 1:
        listy = pseKNCHelper(sequenceX, olinucs, props, Type, k, j, w, geneticMaterial, supInfo)
    else:
        listy = pseKNCHelper(sequenceX, olinucs, props, Type, k, j * gamma, w, geneticMaterial, supInfo)

    return listy


#  getInputSeqs
# ----------------------------------------

#Think about how to return labels and sequences
def getInputSeqs(inputFile):
    label = ''
    labels = []
    sequence = ''
    sequences = []
    tooples = []
    InFileName = inputFile
    InFile = open(InFileName, 'r')
    for line in InFile:
        f = re.search('\>', line)
        g = re.search('(^[A*a*C*c*T*t*G*g*]*$)', line)
        h = re.search('(^[A*a*C*c*G*g*U*u*]*$)', line)
        if f and not (sequence == ''):
            label = line
            labels.append(label)
            sequences.append(sequence)
            sequence = ''
        elif f:
            label = line
            labels.append(label)
        elif g and not (g.group(1) == ''):
            sequence = sequence + g.group(1).upper()
        elif h and not (h.group(1) == ''):
            sequence = sequence + h.group(1).upper()
    sequences.append(sequence)
    for label in labels:
        tooples.append((label, sequences[labels.index(label)]))
    return tooples
    InFile.close()


#  generate_permutations
# ----------------------------------------
# inputs: chars = int
# outputs: all possible oligonucleotides of length k

def generate_permutations(chars):
    allowed_chars = ['A', 'T', 'G', 'C']
    status = []
    for tmp in range(chars):
        status.append(0)
    last_char = len(allowed_chars)
    rows = []
    for x in xrange(last_char ** chars):
        rows.append("")
        for y in range(chars - 1, -1, -1):
            key = status[y]
            rows[x] = allowed_chars[key] + rows[x]
        for pos in range(chars - 1, -1, -1):
            if (status[pos] == last_char - 1):
                status[pos] = 0
            else:
                status[pos] += 1
                break;

    return rows


#  simplePseKNC: Calculates the frequencies of all possible
#  oligonucleotides in the input sequence
# ----------------------------------------
# inputs: inputFile = a string, outputFile = a string, k = int, 
#         formatt = string
# output: nothing, writing to a file

def simplePseKNC(sequenceX, k, formatt, geneticMaterial):
    listy = ''
    # Need to generate all possible oligonucleotides
    olinucz = generate_permutations(k)
    listy = sequenceX
    listz = []
    for olinuc in olinucz:
        freq = listy.count(olinuc)
        freq = int(freq / len(listy) * 1000) / 1000.0
        listz.append(freq)

    return listz

#  main: Gets arguments from the command line, calls either
#  simpePseKNC or pseKNC
# -----------------------------------------
if __name__ == "__main__":
    sequenceX = 'gccggagcgggcttctatttcctcaagggcaacatggtgctgctggacatggctctgcagcggttcgcgatggacctgctggtcaagagaggcttcaccccggtcaccccgccgtacatgatgaaccgcaagtcctacgaaggggtcacagatttagccgacttcgagaaggtcatgtacaagatcgagggcgacgacatgtacctgatcgccacaagcgagcaccccatcggagcaatg'

    Type = 1
    kay = 2
    weight = 0.5
    lam = 1
    formatt = 'csv'
    geneticMaterial = 'DNA'

    props = ['Breslauer_dG', 'Breslauer_dS', 'Electron_interaction', 'Entropy']
    props1 = ['Breslauer_dG', 'Breslauer_dS']
    print(pseKNC(sequenceX.upper(), props1, Type, kay, lam, weight, formatt, geneticMaterial))
    print(simplePseKNC(sequenceX.upper(), 3, formatt, geneticMaterial))
	
