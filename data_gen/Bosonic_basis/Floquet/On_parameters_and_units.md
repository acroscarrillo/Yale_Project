# A word about the parameters and units used for the Floquet code

First things first, as a reminder, the Hamiltonian of interest studied here is 
$$
H(t)/\hbar = \omega_0 a^\dagger a  + g_3(a+a^\dagger)^3 + g_4(a+a^\dagger)^4 \\ - 2i\Omega_1 \cos(\omega_1 t) (a-a^\dagger) - 2i\Omega_2 \cos(\omega_2 t) (a-a^\dagger).
$$
where note that, by dividing by $\hbar$, all the parameters are now in units of frequency $Hz$. Experimentalits like to work with these units, so plots like 
![alt text](experimental_decay_times.png)
are produced. From looking at the plot above, parameters $\epsilon_1$ and $\epsilon_2$ max out at
$$
\frac{\epsilon_1}{2\pi}  = 8MHz,\quad \frac{\epsilon_2}{2\pi}  = 10MHz \quad \text{with} \quad \frac{\omega_0}{2\pi}  = 6GHz.
$$
These parameters are related to the original microscopic parameters via 
$$
\epsilon_1 = 2\Omega_1 \quad \text{and} \quad \epsilon_2 = \frac{2 g_3\Omega_2}{3\omega_0},
$$
where the expression for $\epsilon_2$ was obtained from the paper "for your eyes only" and the expression for $\epsilon_1$ from Rodrigo's mind <span style="color:red">(I need to check these expressions)</span>. Of course $\epsilon_1$ and $\epsilon_2$ are coefficients of an effective (time-independent) Hamiltonian of $H(t)$ above. In many of the papers around the experiment, they work instead with units of $\omega_0$. This is, they work in the dimensionless units
$$
\frac{\epsilon_1}{\omega_0}  = 8/6 \approx 0.0013,\quad \frac{\epsilon_2}{\omega_0} \approx 0.0016 \quad \text{with} \quad \frac{\omega_0}{2\pi}  = 6GHz.
$$
Of course, this redefinition of units translates to $\Omega_1$ and $\Omega_2$ via
$$
\left(\frac{\epsilon_1}{\omega_0} \right)= 2 \left(\frac{\Omega_1}{\omega_0} \right) \quad \text{and} \quad \left(\frac{\epsilon_2}{\omega_0} \right) = \frac{2}{3} \left(\frac{g_3}{\omega_0}\right)\left(\frac{\Omega_2}{\omega_0}\right).
$$
To quote some numbers in these units, in the figure above $\Omega_1$ and $\Omega_2$ max out at
$$
\max\frac{\Omega_1}{\omega_0} = \frac{8/(6\times 10^3)}{2}\approx 0.00066
$$
and
$$
\max \frac{\Omega_2}{\omega_0}  = \frac{3\times10/(6\times10^3)}{2\times0.00075}\approx 3.33,
$$
where $g_3/\omega_0=0.00075$ was an example taken from the paper "for your eyes only".