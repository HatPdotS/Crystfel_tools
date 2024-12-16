from crystfel_tools.crystfel_tools import convert_crystfel_to_mtz


filein = ''                                                     # Input file name
fileout = filein.replace('.hkl','.mtz')                         # Output file name
cell = []                                                       # Unit cell parameters
Symm = ''                                                       # Space group   



convert_crystfel_to_mtz(filein,fileout,cell, Symm)


