__author__ = 'reymolina'

'''
 Ponderations for RNA/DNA and amino acids
'''
import shelve

nitrog_bases = {'A': 0.22,
                'G': 0.24,
                'T': 0.21,
                'C': 0.19,
                'U': 0.21}

amino_acids = {
    'A': 2.07,
    'G': 1.885,
    'L': 1.914,
    'S': 2.09,
    'V': 2.24,
    'T': 2.257,
    'K': 2.254,
    'D': 1.997,
    'I': 2.026,
    'N': 2.071,
    'E': 1.903,
    'P': 1.42,
    'R': 2.211,
    'F': 1.862,
    'Q': 1.86,
    'Y': 2.257,
    'H': 2.041,
    'C': 1.997,
    'M': 1.911,
    'W': 1.886,
    'B': 2.071,
    'Z': 1.86,
    'X': 1.34}

AASeq = 'ALAGLYLEUSERVALTHRLYSASPILEASNGLUPROARGPHEGLNTYRHISCYSMETTRPASXGLXUNK'
AASEQ2 = 'AGLSVTKDINEPRFQYHCMWBZX'
NAcids = 'AGTCU'

def getPond(n, macrotype=0):
    if macrotype == 0:
        k = shelve.open('datapond.db', flag='r')
        #print k['Dname']
        pondname = k['Dname'][n]
        pondvalues = k['Dvalue'][n]
        #print pondname
        #print pondvalues
        k.close()
    else:
        k = shelve.open('dataNCpond.db', flag='r')
        pondname = k['Dname'][n]
        pondvalues = [k['Dvalues'][1]['A'], k['Dvalues'][1]['G'], k['Dvalues'][1]['T'], k['Dvalues'][1]['C'], k['Dvalues'][1]['U']]

        print('Ï was here...')

    return pondname, pondvalues

def getdict (n, macrotype=0):
    ponddict = {}

    if macrotype == 0:
        k = shelve.open('datapond.db', flag='r')
        #print k['Dname']
        pondname = k['Dname'][n]
        pondvalues = k['Dvalue'][n]
        k.close()
        counter = 0
        for i in AASEQ2:
            ponddict[i]= pondvalues[counter]
            counter += 1
    else:
        pondname = 'sum-Amber95'
        ponddict = nitrog_bases


    return pondname, ponddict

if __name__ == "__main__":
    getdict(0, 0)
