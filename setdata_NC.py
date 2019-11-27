__author__ = 'reymolina'

import shelve

NC_Pondx = ['sum-Amber95',
           'first oscillator strength value ( singlet excitation energies)',
           'second oscillator strength value',
            'molar absorption coefficient at 260 nm and PH = 7.0',
           'first (DE1) single excitation energies in eV',
            'second (DE2) single excitation energies in eV']

NC_data = [
    {'A': 0.22, 'G': 0.24, 'U': 0.21, 'T': 0.21, 'C': 0.19 },
    {'A': 0.28, 'G': 0.20, 'U': 0.18, 'T': 0.18, 'C': 0.13},
    {'A': 0.54, 'G': 0.27, 'U': 0.3, 'T': 0.37, 'C': 0.72},
    {'A': 15.4, 'G': 11.7, 'U': 9.9, 'T': 9.2, 'C': 7.5},
    {'A': 4.75, 'G': 4.49, 'U': 4.81, 'T': 4.67, 'C': 4.61},
    {'A': 5.99, 'G': 5.03, 'U': 6.11, 'T': 5.94, 'C': 6.26}]


print ('Ahora voy .....')
s = shelve.open('dataNCpond.db')
try:
    s['Dname'] = NC_Pondx
    s['Dvalues'] = NC_data

finally:
    s.close()


k = shelve.open('dataNCpond.db', flag='r')
print(k['Dname'][1])
print(k['Dvalues'][1]['C'])

k.close()
