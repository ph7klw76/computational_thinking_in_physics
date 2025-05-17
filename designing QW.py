import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e
from scipy.linalg import eigh_tridiagonal

# Physical constants
m_star = 0.066 * m_e    # effective mass inside well
eV_to_J = e
J_to_eV = 1 / e

Delta_Ec_nominal = 0.30  # eV
hv = 0.155               # eV

def ground_state_energy(L_nm, Delta_Ec=Delta_Ec_nominal, dx_nm=0.05):
    """
    Ground state energy (eV) for a finite square well of width L_nm (nm).
    """
    pad_nm = 10.0
    z_nm = np.arange(-pad_nm - L_nm/2, pad_nm + L_nm/2 + dx_nm, dx_nm) # length of quantum well
    N = len(z_nm)
    V = np.where(np.abs(z_nm) <= L_nm/2, 0.0, Delta_Ec)  # eV
    
    dx = dx_nm * 1e-9  # meters
    prefactor = (hbar**2) / (2 * m_star * dx**2) * J_to_eV  # eV
    
    main_diag = 2 * prefactor + V
    off_diag = -prefactor * np.ones(N-1)
    
    eigvals, eigvecs = eigh_tridiagonal(main_diag, off_diag, select='i', select_range=(0,0))
    E1 = eigvals[0]
    psi1 = eigvecs[:,0]
    # normalize
    norm = np.trapz(np.abs(psi1)**2, z_nm*1e-9)
    psi1 /= np.sqrt(norm)
    return E1, z_nm, V, psi1

# Sweep L values
L_values = np.linspace(0.0, 10.0, 1000)
E1_values = np.array([ground_state_energy(L)[0] for L in L_values])

target = Delta_Ec_nominal - hv
tolerance = 0.001  # eV

indices = np.where(np.abs(E1_values - target) <= tolerance)[0]
L_star = L_values[indices[0]] if len(indices)>0 else None

def optimise_L_for_DeltaEc(Delta_Ec_new):
    E1_vals = [ground_state_energy(L, Delta_Ec=Delta_Ec_new)[0] for L in L_values]
    target_new = Delta_Ec_new - hv
    idx = [i for i, E in enumerate(E1_vals) if abs(E - target_new) <= tolerance]
    return (L_values[idx[0]], E1_vals[idx[0]]) if idx else (None, None)

L_minus, E1_minus = optimise_L_for_DeltaEc(0.285)
L_plus,  E1_plus  = optimise_L_for_DeltaEc(0.315)

# Plot E1(L)
plt.figure()
plt.plot(L_values, E1_values)
if L_star is not None:
    plt.scatter([L_star], [E1_values[indices[0]]])
plt.xlabel('Well width L (nm)')
plt.ylabel('Ground state energy E1 (eV)')
plt.title('E1 vs Well Width')
plt.grid(True)
plt.tight_layout()
plt.show()

# Potential and wavefunction for L*
if L_star is not None:
    E1_star, z_nm_star, V_star, psi1_star = ground_state_energy(L_star)
    plt.figure()
    plt.plot(z_nm_star, V_star, label='V(z)')
    scale = Delta_Ec_nominal
    plt.plot(z_nm_star, psi1_star/np.max(np.abs(psi1_star))*scale, label='ψ1 (scaled)')
    plt.axhline(E1_star, linestyle='--', label=f'E1 = {E1_star:.3f} eV')
    plt.xlabel('z (nm)')
    plt.ylabel('Energy / scaled ψ1 (eV)')
    plt.title(f'Potential & Ground State Wavefunction for L* = {L_star:.2f} nm')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

results = {
    "L_star (nm)": L_star,
    "E1(L_star) (eV)": E1_values[indices[0]] if len(indices)>0 else None,
    "L_star (DeltaEc -5%) (nm)": L_minus,
    "L_star (DeltaEc +5%) (nm)": L_plus
}
