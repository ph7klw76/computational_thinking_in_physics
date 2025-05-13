# Designing a GaAs/Al₀.₃Ga₀.₇As Quantum‑Well Infra‑Red Photodetector for 8 µm

![image](https://github.com/user-attachments/assets/14c22695-6238-4439-bb11-2089a26517b8)

Infra‑red QWIPs (quantum‑well infrared photodetectors) detect light by photo‑exciting electrons from the bound ground state in a narrow GaAs quantum well into the continuum of the surrounding AlGaAs barrier, after which an applied bias sweeps them out as photocurrent. Your R & D group wants a single‑well QWIP that peaks at $\lambda = 8.0\mu\text{m}$ (photon energy $h\nu = 0.155\text{eV}$).

## Material Parameters

| Material             | Conduction‑band effective mass $m^*$ | Band gap $E_g$ (300 K) | Conduction‑band offset $\Delta E_c$ |
|----------------------|---------------------------------------|------------------------|--------------------------------------|
| GaAs (well)          | $0.066 m_0$                          | $1.424\text{eV}$     | —                                    |
| Al₀.₃Ga₀.₇As (barrier) | $0.092 m_0$ *(you may assume GaAs mass for a single-mass model)* | $1.925\text{eV}$ | $\Delta E_c \simeq 0.30\text{eV}$ |

Treat the conduction band as a **finite square well**:

$$
V(x) = 
\begin{cases}
0 & \text{if } |x| < \frac{L}{2} \\\\
\Delta E_c & \text{if } |x| \geq \frac{L}{2}
\end{cases}
$$

Assume the material is **n-type**, so only electrons (conduction band) matter.

---

## Tasks

### 1. Write a Python program that, for a given well width $L$ (in metres or nanometres, but be consistent):
- Solves the time‑independent Schrödinger equation for the bound states (use  a finite-difference Hamiltonian method);
- Returns the ground-state energy $E_1(L)$ measured from the bottom of the well.

### 2. Iterate over $L$ to find the smallest width $L^* < 10\,\text{nm}$ for which an 8 µm photon can just promote an electron into the continuum, i.e.

![image](https://github.com/user-attachments/assets/92a3510d-de13-458e-bb16-e4e546254776)


![image](https://github.com/user-attachments/assets/48778b76-2907-4754-977d-780d193fd901)

---

### 3. Visualise
- The potential profile and the normalised ground-state wavefunction for $L^*$
- A graph of $E_1(L)$ versus $L$ from 2 nm to 10 nm, marking the design point $L^*$

---

### 4. Report
- $L^*$ (in nm);
- $E_1(L^*)$ and the residual error (in meV);
- How sensitive the design is to a ±5 % uncertainty in $\Delta E_c$ (re‑optimise for $\Delta E_c = 0.285\text{eV}$ and $0.315\text{eV}$).

## Guide to Solution
## Mathematical Setup

### Effective-Mass Schrödinger Equation

We solve the 1D time-independent effective-mass Schrödinger equation:

$$
\frac{\hbar^2}{2 m^*} \frac{d^2 \psi(z)}{dz^2} + V(z) \psi(z) = E \psi(z)
$$

---

### Potential Profile

A quantum well of width $L$ centered at $z = 0$ with barrier height $\Delta E_c$:

$$
V(z) = 
\begin{cases}
0, & |z| \leq \frac{L}{2} \\\\
\Delta E_c, & |z| > \frac{L}{2}
\end{cases}
$$

---

### Discretization

Define grid points:

$$
z_i = z_0 + i \Delta x, \quad i = 0, \dots, N - 1
$$

Second derivative using central differences:
To learn more calculation of derivatives using central differences method click [here](https://en.wikipedia.org/wiki/Finite_difference)

![image](https://github.com/user-attachments/assets/5052930b-4598-4797-94ed-9dfb927c3a8d)


---

### Hamiltonian Matrix

The Hamiltonian matrix $H$ is tridiagonal with elements:

$$
H_{ii} = \frac{\hbar^2}{m^* (\Delta x)^2} + V(z_i), \quad
H_{i, i \pm 1} = -\frac{\hbar^2}{2 m^* (\Delta x)^2}
$$

---

### Eigenvalue Problem

Solve the matrix equation:

$$
H \Psi = E \Psi
$$

and extract the **lowest eigenvalue** $E_1$ as the ground-state energy.

## Finite-Difference Hamiltonian Matrix

For a grid of $N$ points $z_1, \dots, z_N$ with spacing $\Delta x$, the finite-difference Hamiltonian in the effective-mass approximation is the $N \times N$ tridiagonal matrix:

$$
H =
\begin{pmatrix}
d_1 & -t  & 0   & \cdots & 0 \\
-t  & d_2 & -t  & \ddots & \vdots \\
0   & -t  & d_3 & \ddots & 0 \\
\vdots & \ddots & \ddots & \ddots & -t \\
0 & \cdots & 0 & -t & d_N
\end{pmatrix}
$$

### Matrix Element Definitions

For each grid point $z_i$, define:

- Main diagonal:
  
$$
d_i = H_{ii} = \frac{\hbar^2}{m^* (\Delta x)^2} + V(z_i)
$$

- Off-diagonals:
  
$$
t = -H_{i, i \pm 1} = \frac{\hbar^2}{2 m^* (\Delta x)^2}
$$

So the matrix entries are:

- **Main diagonal** ($i = 1, \dots, N$):

$$
H_{ii} = \frac{\hbar^2}{m^* (\Delta x)^2} + V(z_i)
$$

  where

$$
  V(z_i) =
  \begin{cases}
  0 & \text{inside the well} \\\\
  \Delta E_c & \text{outside the well}
  \end{cases}
$$

- **Off-diagonal elements** ($i = 1, \dots, N - 1$):

$$
H_{i,i+1} = H_{i+1,i} = -\frac{\hbar^2}{2 m^* (\Delta x)^2}
$$

All other entries of $H$ are zero.

---

### Solving the System

Solving the eigenvalue problem

$$
H \psi = E \psi
$$

yields discrete approximations to the **bound-state energies** and **wavefunctions** for your finite square well.

### Libary

You need to use
from [scipy.linalg import eigh_tridiagonal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html).
There are several additional scenarios across quantum mechanics, classical mechanics, and PDE eigenproblems where you can exploit SciPy’s specialized tridiagonal eigensolver such as Quantum harmonic oscillator (1D),Vibrational modes of a linear mass–spring chain, Heat‐equation and etc

An example of part of code is givem below:

```python

dx = dx_nm*1e-9                    # grid spacing in meters
prefactor = (hbar**2)/(2*m_star*dx**2)*J_to_eV  
# here prefactor == a in eV
V = np.where(np.abs(z_nm)<=L/2, 0.0, Delta_Ec)   # potential array V(z_i)
main_diag = 2*prefactor + V        # length-N array [d₁, d₂, …, d_N]
off_diag  = -prefactor * np.ones(N-1)  # length-(N-1) array [–a, –a, …]


from scipy.linalg import eigh_tridiagonal

# select='i', select_range=(0,0) asks only for the lowest eigenpair
eigvals, eigvecs = eigh_tridiagonal(main_diag,
                                     off_diag,
                                     select='i',
                                     select_range=(0,0))

E1   = eigvals[0]       # the ground-state energy
psi1 = eigvecs[:,0]     # the corresponding eigenvector (wavefunction on the grid)
```

Mathematically, you can think of `main_diag` as the diagonal of the Hamiltonian matrix $H$, and `off_diag` as the $\pm 1$ off-diagonals.

Internally, **SciPy** constructs the Hamiltonian as:

![image](https://github.com/user-attachments/assets/51605355-a048-44e6-a08a-c106a0bcc9dc)


This forms a tridiagonal matrix structure:

- `main_diag`: the values along the main diagonal ($H_{ii}$),
- `off_diag`: the values just above and below the diagonal ($H_{i,i\pm1}$).

Then, SciPy applies an **efficient tridiagonal eigensolver**  essentially a specialized version of either the **symmetric QR algorithm** or the **divide-and-conquer algorithm**  to compute the eigenvalues and eigenvectors efficiently.

Selecting Eigenvalues in SciPy

- `select='i'`  
  Tells SciPy to select eigenvalues **by index** rather than computing all of them or selecting by value.

- `select_range=(0, 0)`  
  Together with `select='i'`, this specifies the **inclusive index range** of the eigenvalues to return.  
  Since SciPy uses **0-based indexing**, `(0, 0)` means:

“Return only the eigenvalue whose index is 0,”  
i.e., the **smallest (ground-state) eigenvalue** $E_1$ and its corresponding **eigenvector** $\psi_1$.

You can learn more of the above by clicking [here](https://medium.com/modern-physics/finite-difference-solution-of-the-schrodinger-equation-c49039d161a8)

What you should expect

![image](https://github.com/user-attachments/assets/d75496b6-8db7-4e53-9f48-98887795148a)

