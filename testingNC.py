import shelve

k = shelve.open('dataNCpond.db', flag='r')
print(k['Dname'])
b = k['Dvalues'][0]
print('adenine  {}'.format(b['A']))
print(k['Dvalues'][0])


