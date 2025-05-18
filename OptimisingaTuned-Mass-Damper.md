# Optimising a Tuned-Mass-Damper (TMD) for a Pedestrian Suspension Bridge Deck

When several people walk in step, the horizontal forcing is nearly sinusoidal at $\approx 1.6\,\text{–}\,2.2\ \text{Hz}$ (the cadence of human walking). If the bridge’s first lateral natural frequency sits in this band, resonance can produce deck accelerations that feel unsafe or even cause panic [Millennium Bridge,London,2000](https://www.youtube.com/watch?v=g37pKBl3DfE).

A a smaller “parasitic” mass connected to the deck by a spring and dashpot that absorbs energy near resonance and returns it as heat. Choosing its parameters (mass ratio, stiffness, damping) is a textbook application of forced, damped, coupled oscillations and normal-mode theory.

Light steel footbridges often have a first lateral natural frequency close to the cadence of walking crowds ($\approx 1.6\,\text{–}\,2.2\ \text{Hz}$). Synchronised pedestrians can excite a sinusoidal horizontal force that makes the deck vibrate uncomfortably. A Tuned-Mass Damper (TMD), a small secondary mass connected by a spring-dashpot, re-distributes the vibratory energy and keeps accelerations within the ISO 10137 lateral comfort limit of:
$0.2\ \text{m}\ \text{s}^{-2}$ (RMS)

RMS here stands for root-mean-square, a statistical measure of the magnitude of a varying quantity. In our vibration context, we use it to quantify the “average” acceleration level of the deck over a period of time, regardless of sign (upward or downward).

Concretely, if $a(t)$ is the instantaneous vertical acceleration of the deck, then over a time window $[t_1,\, t_2]$ the RMS acceleration is defined as:

$$
a_{\text{RMS}} = \sqrt{ \frac{1}{t_2 - t_1} \int_{t_1}^{t_2} [a(t)]^2 \, dt }
$$

Unlike a simple arithmetic mean, RMS accounts for both positive and negative swings by squaring, so it measures the energy content of the vibration.

It gives a single scalar “equivalent” acceleration level that you can compare directly against comfort or serviceability thresholds.

# Vibration Model of the Bridge

## Single-Degree-of-Freedom Approximation

Before adding a TMD, we model the bridge’s fundamental vertical mode as a single-degree-of-freedom (SDOF) system:

![image](https://github.com/user-attachments/assets/d0d85ba0-6094-406b-ad8d-d66e2787ce3a)

$$
m\ddot{y}_1 + c\dot{y}_1 + ky_1 = F(t)
$$

- $y_1(t)$: vertical displacement of the bridge (positive downward)  
- $m$: modal mass (kg)  
- $k$: modal stiffness (N/m)  
- $c$: inherent damping (N·s/m), with damping ratio

$$
\zeta = \frac{c}{2m\omega_n}
$$

- $F(t)$: external force (pedestrians, wind), positive downward  
- $\omega_n = \sqrt{\frac{k}{m}}$: natural frequency (rad/s)

Equation (1) follows directly from summing forces on mass $m$: inertial force $m\ddot{y}_1$, damping $c\dot{y}_1$, and spring $ky_1$, balanced by $F(t)$.


## Two-Degree-of-Freedom Model with TMD

We introduce a secondary mass $m_d$ (the TMD) that can move vertically relative to the ground. Define:

- $y_1(t)$: vertical displacement of the primary (bridge) mass $m$  
- $y_2(t)$: vertical displacement of the secondary (TMD) mass $m_d$

Both displacements are measured positive downward from their static equilibrium positions.

![image](https://github.com/user-attachments/assets/a05b043a-c869-40ca-a44b-969a3b2f1006)

For mass $m$, the forces in the vertical direction are:

- Spring $k$ to ground: $-k y_1$  
- Damper $c$ to ground: $-c \dot{y}_1$  
- TMD spring $k_d$ between $m$ and $m_d$: force on $m$ is $-k_d (y_1 - y_2)$  
- TMD damper $c_d$ between $m$ and $m_d$: force on $m$ is $-c_d (\dot{y}_1 - \dot{y}_2)$  
- External force: $+F(t)$

Summing these and equating to $m \ddot{y}_1$ gives:

$$
m\ddot{y}_1 + c\dot{y}_1 + ky_1 + c_d(\dot{y}_1 - \dot{y}_2) + k_d(y_1 - y_2) = F(t)
$$

For mass $m_d$, the only vertical forces are from the TMD elements:

- Spring $k_d$: $-k_d (y_2 - y_1)$  
- Damper $c_d$: $-c_d (\dot{y}_2 - \dot{y}_1)$

No direct external excitation acts on $m_d$, so:

$$
m_d\ddot{y}_2 + c_d(\dot{y}_2 - \dot{y}_1) + k_d(y_2 - y_1) = 0                       
$$


example of tuned-mass damper under a bridge

![image](https://github.com/user-attachments/assets/aa4f75b9-a018-4e49-bcaf-525128290555)


To understand more about tuned-mass damper click [here](https://www.youtube.com/watch?v=vLaFAKnaRJU)
# Compact Matrix Form

We collectthe last 2 equations  into a matrix equation. Let

![image](https://github.com/user-attachments/assets/bfd6bf71-2a01-48e2-9ca3-120dcfd79fae)


and external force vector ![image](https://github.com/user-attachments/assets/d53dad4b-528f-4de2-bff0-0ecea25e2507)
Then,

![image](https://github.com/user-attachments/assets/d22e484e-1a04-46a7-8bf4-713f43982629)

This succinct form highlights the coupling between the two masses via off-diagonal terms in $C$ and $K$.

# Harmonic Response and Transfer Function

Assume sinusoidal excitation $F(t) = F_0 e^{i\omega t}$ and seek steady-state $y(t) = Y e^{i\omega t}$. Substitution into (4) yields:

![image](https://github.com/user-attachments/assets/fa802e67-3b30-48bd-8ee6-7b63b3a51109)


The frequency-response function from force to primary displacement is:

![image](https://github.com/user-attachments/assets/0c966b3c-03f2-41b4-b473-e38401ff5258)

Where $(\cdot)_{1,1}$ denotes the top-left entry of the inverse matrix. Plotting $|H(\omega)|$ reveals how the TMD splits and lowers the resonance peak.

# Optimum Tuning: Den Hartog’s Equal-Peak Method

Define the primary natural frequency

$$
\omega_n = \sqrt{\frac{k}{m}}
$$

and the TMD’s uncoupled frequency

$$
\omega_d = \sqrt{\frac{k_d}{m_d}}.
$$

Let the mass ratio be

$$
\mu = \frac{m_d}{m} \quad (\mu \ll 1\ \text{typical}).
$$

If we choose

$$
r = \frac{\omega_d}{\omega_n}, \quad
\zeta_d = \frac{c_d}{2 m_d \omega_d},
$$

then Den Hartog showed that the equal-peak (minimum-peak) tuning occurs at:

$$
r_{\text{opt}} = \frac{1}{1 + \mu}, \quad
\zeta_{d,\text{opt}} = \frac{3\mu}{8(1 + \mu)}. 
$$

With these, the two resonance peaks in $|H(\omega)|$ are equal and as small as possible.

