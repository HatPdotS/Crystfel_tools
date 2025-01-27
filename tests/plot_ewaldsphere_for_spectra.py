import h5py
from crystfel_tools.handling import fast_math
import numpy as np
from matplotlib import pyplot as plt
keyy = 'SARFE10-PSSS059:SPECTRUM_Y'
keyx = 'SARFE10-PSSS059:SPECTRUM_X'
with h5py.File('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/run_001192.BSREAD.h5', 'r') as f:
    datax = f['/data/'+keyx]['data'][0]
    datay = f['/data/'+keyy]['data'][0]
    datay = datay - datay[:100].mean()


cell = [61.51,91.01,151.11,90,90,90]
def convert_energy_to_ech(Energy):
    c = 3e8 * 1e10
    h = 6.62607015e-34
    e = Energy * 1.60217662e-19
    return e / (h * c)
size = 1000
x = np.linspace(0,1,size)
y = 0
z = np.linspace(1.9,2,size)
grid = np.meshgrid(x,y,z)
grid = np.array(grid).reshape(1,3,-1)
res = np.zeros(size*size)
spectral_stepsize = convert_energy_to_ech((datax[1:] - datax[:-1]).mean())
spectra_ech = convert_energy_to_ech(datax)
for i in range(datax.shape[0]):
    sol = fast_math.get_delta(grid,0,0,convert_energy_to_ech(datax[i]))

    res[sol.flatten()<spectral_stepsize*10] += datay[i]

plt.imshow(res.reshape(size,size))   

plt.savefig('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/test.png')
plt.close()
def calc_diff(hkl0,ech):
    diff = ((hkl0[:,0,:]**2 + hkl0[:,1,:]**2 + hkl0[:,2,:]**2)/(2*hkl0[:,2,:])) - ech
    return  diff
sol = calc_diff(grid,convert_energy_to_ech(spectra_ech.mean()))
idx = np.argmin(np.abs(sol.reshape(-1,1) - (spectra_ech.reshape(1,-1) - spectra_ech.mean())),axis=1)
res += datay[idx]
plt.imshow(res.reshape(size,size))
plt.savefig('/das/work/p17/p17490/Peter/Library/crystfel_tools/test_data/test2.png')
plt.close()
