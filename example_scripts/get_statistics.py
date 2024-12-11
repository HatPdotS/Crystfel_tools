import crystfel_tools.crystfel_tools as cft

cell_path = ''  # set the path to the cell file
hkl_path = ''   # set the path to the hkl file
pg = ''         # set the pointgroup of your crystals
resmax = 2      # maximum resolution in Angstrom (low)
resmin = 10     # minimum resolution in Angstrom (high)
nshells = 20    # number of shells



cft.get_statistics(hkl_path,pg,cell_path,resmax=resmax,resmin=resmin,nshells=nshells)