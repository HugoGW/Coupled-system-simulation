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
    \begin{array}{ll}
        \displaystyle \frac{dn_1}{dt} = F_A - F_D - 2D_1n_1^2 - D_1n_1n - 2\alpha F_A n_1 - \alpha F_A n \\
        \displaystyle \frac{dn}{dt} = D_1n_1^2 + \alpha F_A n_1
    \end{array}
$$


where $n_1$ is the amount of monomers and $n$ is the amount of island, this system is called the rate equation model. These equations are solved with the $\texttt{odeint}$ integrator for the $\texttt{scipy}$ library

We now consider a grid of size (100×100) without edge. In other words, it has the same topological characteristics as the PACMAN game, a monomer that exits to the right reappears to the left or, if it disappears at the bottom, reappears at the top. The grid can be seen as a succession of grids of the same type:

![image](https://github.com/user-attachments/assets/595314fa-83c2-4296-8da7-cf4d0c4a9c69)

From this image, it is clear that the monomers reappear on the other side when they leave the grid.
This grid configuration allows the monomers to be conserved, so there is no loss, and it establishes a link between the monomers located at the ends of the grid, thus encouraging their interactions

In the code, we have the parameters that initiate the simulation.

        Fa = 1             # Deposition rate (new monomer per step)
        Ndif = 5           # Number of diffusion steps for each monomer
        time_steps = 2100  # Total number of time steps
        grid = np.zeros((L, L))  # Initialize an empty grid (0 = empty, 1 = monomer, 2 = island)

        # Set parameters for differential equations
        Fd = 0  # Example value for disappearance rate (death rate) of monomers
        D1 = 0.00005  # Example value for D1 coefficient (aggregation rate)
        a = 0.0085  # Example value for a coefficient (growth rate of islands)

And we obtain coherent results from random diffusion :

![image](https://github.com/user-attachments/assets/3c9f7dfa-52fb-4f7c-802a-242e7a5d494f)






