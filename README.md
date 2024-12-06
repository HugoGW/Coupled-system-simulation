# Coupled-system-simulation
I've made a simulation and a fit of the dynamic of a coupled equation (involving monomers and islands)

We consider a 2D system made of monomers (that will be represented by dots later for simulations) and "islands" that are aggrations of monomers (made in different colors to differentiate them) :
![image](https://github.com/user-attachments/assets/b2284e51-5caf-42cc-8b19-e06a9dbed980)

Monomers can freely move, their motions will be considred as a random diffusion and they're appearing once per second. Where 2 monomers collide or are next to each other, they agglomerate and become a motionless island.
![image](https://github.com/user-attachments/assets/13f73e16-c6e2-4b1d-b576-c3bc7d250dfb)

When a monomer meet an island, it becomes part of the island and make it bigger.
![image](https://github.com/user-attachments/assets/a3996be7-6fa4-47e8-a38c-09da4d5ffa04)

It exists a system of coupled equations that allows us to know the amount of monomers and island according to time :

$$
\left\{
    \begin{array}{ll}
        \displaystyle \frac{dn_1}{dt} = F_A - F_D - 2D_1n_1^2 - D_1n_1n - 2\alpha F_A n_1 - \alpha F_A n ~~~~ (1)\\
        \displaystyle \frac{dn}{dt} = D_1n_1^2 + \alpha F_A n_1 ~~~~ (2)
    \end{array}
\right.
$$

where $n_1$ is the amount of monomers and $n$ is the amount of island, this system is called the rate equation model. These equations are solved with the \texttt{odeint} integrator for the \texttt{scipy} library





