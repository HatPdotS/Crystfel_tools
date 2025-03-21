import numba
import numpy as np



def calc_cc(set1,set2):
    """
    Calculate the correlation coefficient between two sets of data
    """
    return np.corrcoef(set1.flatten(),set2.flatten())[0,1]


def calc_rsplit_weighted(set1,set2,weights=None):
    """
    Calculate the Rsplit between two sets of data
    """
    if weights is None:
        weights = np.ones_like(set1)
    set1 = set1 / set1.mean()
    set2 = set2 / set2.mean()
    return np.sqrt(np.sum((set1 - set2)**2 * weights)) / weights.sum()


def calc_mean(data,weights):
    sum_weigths = np.sum(weights,axis=1).todense()
    sum_weigths[sum_weigths == 0] = 1
    return np.sum(data,axis=1).todense() / sum_weigths

def write_out_half_datasets(data, weigths):
    half1 = data[:,:data.shape[1]//2]

    half2 = data[:,data.shape[1]//2:]

    weigths1 = weigths[:,:weigths.shape[1]//2]
    weigths2 = weigths[:,weigths.shape[1]//2:]
    half1 = calc_mean(half1,weigths1)
    half2 = calc_mean(half2,weigths2)
    return half1,half2

def get_delta(lattice,alpha,beta,ech):
    return np.abs(np.sqrt((lattice[:,0,:] - ech *np.sin(alpha) * np.cos(beta))**2 + (lattice[:,1,:] - ech * np.sin(alpha) * np.sin(beta))**2 + (lattice[:,2,:] - ech * np.cos(alpha))**2) - ech)

def energy_at_S2(lattice,alpha,beta,ech,deltaech):
    delta = get_delta(lattice,alpha,beta,ech)
    energy_at = np.exp(-(delta**2)/deltaech**2)
    return energy_at

def calc_rec_space_vector(cell):
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alpha = np.deg2rad(cell[3])
    beta = np.deg2rad(cell[4])
    gamma = np.deg2rad(cell[5])
    va = np.array([a,0,0])
    vb = np.array([b * np.cos(gamma),b * np.sin(gamma),0])
    Vc = np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))
    vc = np.array([c * np.cos(beta),c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),c * Vc])
    V = np.dot(va,np.cross(vb,vc))
    astar = np.cross(vb,vc) / V
    bstar = np.cross(vc,va) / V
    cstar = np.cross(va,vb) / V
    return np.vstack((astar,bstar,cstar)).T

def log_energy_at_S2(lattice,alpha,beta,ech,deltaech):
    delta = get_delta(lattice,alpha,beta,ech)
    return -(delta**2)/deltaech**2

def rotate_reciprocal_space(hkl,alpha,beta):
    rotation_matrix = np.array([[np.cos(alpha) * np.cos(beta), -np.sin(beta), np.sin(alpha) * np.cos(beta)],
                     [np.cos(alpha) * np.sin(beta), np.cos(beta), np.sin(alpha) * np.sin(beta)],
                     [-np.sin(alpha), 0, np.cos(alpha)]])
    h = hkl[:,0] * rotation_matrix[0,0] + hkl[:,1] * rotation_matrix[1,0] + hkl[:,2] * rotation_matrix[2,0]
    k = hkl[:,0] * rotation_matrix[0,1] + hkl[:,1] * rotation_matrix[1,1] + hkl[:,2] * rotation_matrix[2,1]
    l = hkl[:,0] * rotation_matrix[0,2] + hkl[:,1] * rotation_matrix[1,2] + hkl[:,2] * rotation_matrix[2,2]
    hkl[:,0] = h
    hkl[:,1] = k
    hkl[:,2] = l
    return hkl

def calc_energy_S2_analytical(hkl0,alpha,beta,ech,deltaech):
    hkl0 = rotate_reciprocal_space(hkl0,alpha,beta)
    hkl0[hkl0 == 0] = 1e-10
    diff = ((hkl0[:,0,:]**2 + hkl0[:,1,:]**2 + hkl0[:,2,:]**2)/(2*hkl0[:,2,:])) - ech
    return np.exp(-diff**2 / deltaech**2) / np.abs(hkl0[:,2,:])

def calc_energy_S2_analytical_log(hkl0,alpha,beta,ech,deltaech):
    hkl0 = rotate_reciprocal_space(hkl0,alpha,beta)
    hkl0[hkl0 == 0] = 1e-10
    diff = ((hkl0[:,0,:]**2 + hkl0[:,1,:]**2 + hkl0[:,2,:]**2)/(2*hkl0[:,2,:])) - ech
    return -diff**2 / deltaech**2 - np.log(np.abs(hkl0[:,2,:]))

def make_sampling_box(hkl,cell,mesh_size = 0.2,mesh_nr=11):
    rec = calc_rec_space_vector(cell)
    mesh_size = np.array((0.2,0.2,0.2))
    lattice_points = convert_to_rec_space(hkl,rec)
    mesh_size = convert_to_rec_space(mesh_size,rec).reshape(3)
    # 4 dimensions 1 different crystals 2 h,k,l 3 values values start with
    hkl0s = []
    for triplet in lattice_points:
        h0 = np.linspace(triplet[0]-mesh_size[0],triplet[0]+mesh_size[0],mesh_nr)
        k0 = np.linspace(triplet[1]-mesh_size[1],triplet[1]+mesh_size[1],mesh_nr)
        l0 = np.linspace(triplet[2]-mesh_size[2],triplet[2]+mesh_size[2],mesh_nr)
        grid = np.meshgrid(h0,k0,l0)
        hkl0s.append(np.array(grid).reshape(1,3,-1))
    return np.concatenate(hkl0s,axis=0), lattice_points

def convert_to_rec_space(hkl,rec_vectors):
    hkl = hkl.reshape(-1,3,1)
    rec_vectors = rec_vectors.reshape(1,3,3)
    return np.sum(hkl * rec_vectors,axis=1)

def intensity_of_all_non_observed(hkl_not_observed,alpha,beta,ech,deltaech,cell):
    rec = calc_rec_space_vector(cell)
    lattice_points = convert_to_rec_space(hkl_not_observed,rec).reshape(-1,3,1)
    return np.sum(energy_at_S2(lattice_points,alpha,beta,ech,deltaech))

def I_at_hkl0(hkl0,hkl,Bi):
    sigma = np.sum(hkl**2 * Bi,axis=1)
    return np.exp(-np.sum((hkl0-hkl)**2,axis=-2) / sigma) / sigma


def calc_rho(hkl,alpha,beta,ech,deltaech,Bi,cell):
    lattice, lattice_points = make_sampling_box(hkl,cell)
    rho = np.sum(I_at_hkl0(lattice,lattice_points.reshape(-1,3,1),Bi.reshape(1,3,1)) * energy_at_S2(lattice,alpha,beta,ech,deltaech),axis=-1)
    rho/= np.mean(rho)
    return rho


def plot_I_at_hkl0(hkl,Bi,cell):
    lattice, lattice_points = make_sampling_box(hkl,cell)
    return I_at_hkl0(lattice,lattice_points.reshape(-1,3,1),Bi.reshape(1,3,1))

def calc_log_xsphere(hkl,alpha,beta,ech,deltaech,cell):
    _, lattice_points = make_sampling_box(hkl,cell)
    return calc_energy_S2_analytical_log(lattice_points.reshape(-1,3,1),alpha,beta,ech,deltaech).reshape(hkl.shape[0])

def residual_xpshere(factors,hkl,I,I_ref,cell,ech,deltaech):
    res = np.log(I) - np.log(I_ref) - calc_log_xsphere(hkl,*factors[[0,1]],ech,deltaech,cell) - factors[0] + np.sum(hkl**2,axis=1) * factors[3]
    return np.sum(res ** 2)

def residual_simple(factors,hkl,I,I_ref,cell):
    rec = calc_rec_space_vector(cell)
    lattice = convert_to_rec_space(hkl,rec)
    res = I*(np.log(I) - (np.log(I_ref) + np.log(factors[0]) - np.sum(lattice**2,axis=1) * factors[1]))
    return np.sum(res ** 2)

def apply_scaling_simple(factors,hkl,I,cell):
    rec = calc_rec_space_vector(cell)
    lattice = convert_to_rec_space(hkl,rec)
    p = np.exp(-np.sum(lattice**2,axis=1) * factors[1]) * np.exp(factors[0])
    p = p/p.mean()
    mask = p > 0.1
    I_scaled = I[mask] / p[mask]
    return I_scaled, p[mask], mask

def apply_scaling_xsphere(factors,hkl,I,cell,ech,deltaech):
    p = np.exp(calc_log_xsphere(hkl,*factors[[1,2]],ech,deltaech,cell))
    p = p/p.mean()
    mask = p > 0.1 
    Bterm = np.exp(np.sum(hkl**2,axis=1) * factors[3])[mask]
    I_scaled = I[mask] / np.exp(factors[0]) / p[mask] / Bterm
    weights = factors[0] * p[mask] * Bterm
    weights = weights/weights.mean()
    return I_scaled, weights, mask

def calc_logrho(hkl,alpha,beta,ech,deltaech,Bi,cell):
    rho = calc_rho(hkl,alpha,beta,ech,deltaech,Bi,cell)
    msk = rho == 0
    rho[msk] = 1
    logrho = np.log(rho)
    logrho[msk] = 1000
    return logrho

def calc_rho_new(hkl,alpha,beta,ech,deltaech,Bi,cell):
    lattice, lattice_points = make_sampling_box(hkl,cell)
    rho = np.sum(I_at_hkl0(lattice,lattice_points.reshape(-1,3,1),Bi) * calc_energy_S2_analytical(lattice,alpha,beta,ech,deltaech),axis=-1)
    rho/= np.max(rho)
    return rho

def calc_logrho_new(hkl,alpha,beta,ech,deltaech,Bi,cell):
    rho = calc_rho_new(hkl,alpha,beta,ech,deltaech,Bi,cell)
    msk = rho == 0
    rho[msk] = 1
    logrho = np.log(rho)
    logrho[msk] = -1000
    return logrho


def calc_residual_new(factors,hkl,I,I_ref,cell,ech,deltaech):
    Bi = factors[3]
    weights = I/I.mean()
    return np.sum(weights * (np.log(I) - (np.log(I_ref) + factors[0] + calc_logrho_new(hkl,*factors[[1,2]],ech,deltaech,Bi,cell))) ** 2)

def residual_with_rho(factors,hkl,hkl_not_observed,I,I_ref,cell):
    Bi = factors[5:8]
    return np.sum((np.log(I) - np.log(I_ref) - np.log(factors[0]) - calc_logrho(hkl,*factors[[1,2,3,4]],Bi,cell)) ** 2 ) + intensity_of_all_non_observed(hkl_not_observed,factors[1],factors[2],ech,deltaech,cell)

def residual_with_rho_isotropic(factors,hkl,hkl_not_observed,I,I_ref,cell,ech,deltaech):
    Bi = np.array((factors[3],factors[3],factors[3]))
    return np.sum((np.log(I) - (np.log(I_ref) + factors[0] + calc_logrho(hkl,*factors[[1,2]],ech,deltaech,Bi,cell))) ** 2)

def evaluate_crystal(factors,hkl,I,ech,delatech,cell,Bi=None):
    if Bi is None:
        Bi = factors[5:8]
    rho = calc_rho(hkl,*factors[[1,2]],ech,delatech,Bi,cell)
    mask = rho > 0.1
    I = I[mask]
    rho = rho[mask]
    weights = factors[0] * rho
    weights = weights/weights.mean()
    I_scaled = I / factors[0] / rho
    return  I_scaled, weights, mask

def evaluate_crystal_new(factors,hkl,I,cell,ech,deltaech,rho_cut=0.1):
    Bi = factors[3]
    rho = calc_rho_new(hkl,*factors[[1,2]],ech,deltaech,Bi,cell)
    mask = rho > rho_cut
    I = I[mask]
    rho = rho[mask]
    weights = factors[0] * rho
    weights = weights/weights.mean()
    I_scaled = I / factors[0] / rho
    return  I_scaled, weights, mask

def evaluate_crystal_xsphere(factors,hkl,I,cell):
    return I / factors[0] / np.exp(calc_log_xsphere(hkl,*factors[[1,2,3,4]],cell))

def reciprocal_basis_matrix(unit_cell):
    # Extract unit cell parameters
    a, b, c, alpha, beta, gamma = unit_cell
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    # Compute real-space basis vectors
    cos_alpha, cos_beta, cos_gamma = np.cos(alpha), np.cos(beta), np.cos(gamma)
    sin_gamma = np.sin(gamma)
    volume = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)
    a_vec = np.array([a, 0, 0])
    b_vec = np.array([b * cos_gamma, b * sin_gamma, 0])
    c_vec = np.array([
        c * cos_beta,
        c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
        c * volume / sin_gamma
    ])
    # Compute reciprocal basis vectors
    volume_real = np.dot(a_vec, np.cross(b_vec, c_vec))
    a_star = np.cross(b_vec, c_vec) / volume_real
    b_star = np.cross(c_vec, a_vec) / volume_real
    c_star = np.cross(a_vec, b_vec) / volume_real
    # Assemble reciprocal basis matrix
    return np.array([a_star, b_star, c_star])

def get_scattering_vectors(hkl, unit_cell):
    recB = reciprocal_basis_matrix(unit_cell)
    hkl = np.array(hkl)  # Ensure hkl is a numpy array
    s = np.dot(hkl,recB)
    return s

def get_resolution(hkl, unit_cell):
    s = get_scattering_vectors(hkl, unit_cell)
    return 1 / np.sum(s**2, axis=1)**0.5