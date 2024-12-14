import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import label
from scipy.integrate import odeint
from scipy.optimize import curve_fit

# make sure you installed all the libraries


# Define model parameters
L = 100            # Grid size (L x L)
Fa = 1             # Deposition rate (new monomer per step)
Ndif = 5           # Number of diffusion steps for each monomer
time_steps = 2100  # Total number of time steps
grid = np.zeros((L, L))  # Initialize an empty grid (0 = empty, 1 = monomer, 2 = island)

# Store data for plotting
monomer_counts = []  # List to store the number of monomers at each time step
independent_island_counts = []  # List to store the number of independent islands at each time step
island_percentages = []  # List to store the percentage of grid covered by islands

# Set parameters for differential equations
Fd = 0  # Example value for disappearance rate (death rate) of monomers
D1 = 0.00005  # Example value for D1 coefficient (aggregation rate)
a = 0.0085  # Example value for a coefficient (growth rate of islands)

# Define the system of differential equations for monomer and island dynamics
def system(y, t, Fa, Fd, D1, a):
    n1, n = y  # n1: number of monomers, n: number of islands
    # Differential equations for monomer and island evolution
    dn1_dt = Fa - Fd - 2*D1*n1**2 - D1*n1*n - 2*a*Fa*n1 - a*Fa*n
    dn_dt = D1*n1**2 + a*Fa*n1
    return [dn1_dt, dn_dt]

# Initial conditions for the differential equations
n1_0 = 0  # Initial value for n1 (number of monomers) at time t=0
n_0 = 0   # Initial value for n (number of islands) at time t=0
y0 = [n1_0, n_0]  # Initial state vector

# Time points for theoretical solution of the differential equations
t_theoretical = np.arange(0, time_steps * 30, 0.1)  # Time points for evaluation

# Solve the system of differential equations to obtain the theoretical curve
solution_theoretical = odeint(system, y0, t_theoretical, args=(Fa, Fd, D1, a))

# Function to deposit a monomer at a random empty site on the grid
def deposit_monomer(grid):
    empty_cells = np.argwhere(grid == 0)  # Find all empty cells on the grid
    if len(empty_cells) > 0:  # If there are empty cells available
        new_monomer_pos = empty_cells[np.random.choice(len(empty_cells))]  # Select a random empty cell
        grid[new_monomer_pos[0], new_monomer_pos[1]] = 1  # Place a monomer at the selected position

# Function to diffuse monomers randomly
def diffuse_monomer(grid, pos):
    L = grid.shape[0]  # Get the size of the grid
    x, y = pos  # Current position of the monomer
    # List of possible movements: up, down, left, right
    moves = [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]
    moves = [(i % L, j % L) for i, j in moves]  # Apply periodic boundary conditions at grid edges
    np.random.shuffle(moves)  # Randomize the order of movements
    for new_pos in moves:  # Try each possible movement
        if grid[new_pos[0], new_pos[1]] == 0:  # Move only if the new site is empty
            grid[x, y] = 0  # Clear the current position
            grid[new_pos[0], new_pos[1]] = 1  # Move the monomer to the new position
            return new_pos  # Return the new position
    return pos  # If no movement is possible, stay at the current position

# Function to check if aggregation should occur (monomer joins an island)
def check_aggregation(grid, pos):
    L = grid.shape[0]  # Get the size of the grid
    x, y = pos  # Current position of the monomer
    # List of adjacent sites: up, down, left, right
    moves = [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]
    moves = [(i % L, j % L) for i, j in moves]  # Apply periodic boundary conditions
    for new_pos in moves:
        if grid[new_pos[0], new_pos[1]] == 1 or grid[new_pos[0], new_pos[1]] == 2:  # Check if neighboring site is a monomer or island
            grid[x, y] = 2  # Convert the monomer to an island
            return True  # Aggregation occurred
    return False  # No aggregation

# Function to count the number of independent islands on the grid
def count_independent_islands(grid):
    labeled_array, num_features = label(grid == 2)  # Label connected components (islands)
    return num_features  # Return the number of independent islands

# Function to calculate the total area covered by islands
def calculate_island_area(grid):
    return np.sum(grid == 2)  # Count the number of cells occupied by islands

# Lists to store data for fitting
points_x = []  # Store time steps
points_y_n1 = []  # Store monomer counts for fitting
points_y_n = []  # Store island counts for fitting

# Function to update the simulation at each frame
def update(frame):
    global grid, popt_n1  # Make grid and fitting parameters accessible
    
    # Calculate the total area covered by islands and check condition
    island_area = calculate_island_area(grid)  # Calculate the current island area
    total_area = L * L  # Total grid area
    island_percentage = (island_area / total_area) * 100  # Percentage of grid covered by islands
    island_percentages.append(island_percentage)  # Store the percentage for plotting
    
    if island_percentage < 20:  # If less than 20% of the area is occupied by islands
        # Step 1: Deposit new monomers
        for _ in range(Fa):  # Loop to deposit the number of monomers equal to the deposition rate
            deposit_monomer(grid)  # Call the function to deposit a monomer at a random empty site

    # Step 2: Diffuse the monomers
    monomers = np.argwhere(grid == 1)  # Get all positions of the monomers in the grid
    for monomer in monomers:  # Loop over each monomer's position
        pos = tuple(monomer)  # Convert the position array to a tuple for easy access
        for _ in range(Ndif):  # Loop for the number of diffusion steps defined
            if check_aggregation(grid, pos):  # Check if aggregation occurs at the current position
                break  # If aggregation occurs, stop diffusion for this monomer
            pos = diffuse_monomer(grid, pos)  # Diffuse the monomer to a new position

    # Count the number of monomers and independent islands at this time step
    n1 = np.sum(grid == 1)  # Count the number of monomers present in the grid
    num_islands = count_independent_islands(grid)  # Count the number of independent islands in the grid
    
    # Save the counts for plotting
    monomer_counts.append(n1)  # Append the current number of monomers to the list
    independent_island_counts.append(num_islands)  # Append the current number of islands to the list

    # Update points for curve fitting
    points_x.append(frame)  # Append the current frame number to the x points for fitting
    points_y_n1.append(monomer_counts[frame])  # Append the number of monomers for the current frame
    points_y_n.append(independent_island_counts[frame])  # Append the number of islands for the current frame

    # Fit every 100 steps only
    if frame % 2000 == 0 and len(points_x) >= 5:  # Fit only if there are enough points and it's a multiple of 100
        p0 = [D1, a]  # Initial guess for D1 and a parameters for the fit
        try:
            # Perform curve fitting to adjust D1 and a
            popt_n1, _ = curve_fit(lambda t, D1, a: odeint(system, y0, t, args=(Fa, Fd, D1, a))[:, 0],
                                   points_x, points_y_n1, p0=p0)
            
            # Generate the fitted theoretical curve for n1
            y_fit_n1 = odeint(system, y0, t_theoretical, args=(Fa, Fd, popt_n1[0], popt_n1[1]))[:, 0]
            
            # Update the fitting line for n1 (dashed line)
            fit_line_n1.set_data(t_theoretical, y_fit_n1)
            
            # Generate the fitted theoretical curve for n
            y_fit_n = odeint(system, y0, t_theoretical, args=(Fa, Fd, popt_n1[0], popt_n1[1]))[:, 1]
            
            # Update the fitting line for n (dashed line)
            fit_line_n.set_data(t_theoretical, y_fit_n)
            
            # Update the text displays for D1 and a after each fitting
            d1_text.set_text(f"$D_1$: {popt_n1[0]:.5f}")  # Display the fitted D1 value
            a_text.set_text(f"$\\alpha$: {popt_n1[1]:.5f}")  # Display the fitted a value
            
        except RuntimeError:
            pass  # In case of fitting error, do nothing

    # Update the visualization
    mat.set_data(grid)  # Update the grid visualization for the current frame
    return [mat]  # Return the updated grid for the animation

# Configure the visualization
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))  # Create a 1x3 grid for subplots
mat = ax1.matshow(grid, cmap='coolwarm', vmin=0, vmax=2)  # Show the grid as a heatmap
ax1.set_title("Dynamics of Monomers and Islands", fontsize=10)  # Set title for the grid plot
ax1.set_xticks([])  # Hide x ticks
ax1.set_yticks([])  # Hide y ticks

# Plot for the number of monomers
line1, = ax2.plot([], [], color='royalblue', label='Simulated Monomers', alpha=0.9)  # Line for simulated monomers
ax2.plot(1, -100, '--k', label='Theoretical Fit', lw=2)  # Placeholder for theoretical fit line
monomer_text = ax2.text(0.7, 0.77, '', transform=ax2.transAxes, fontsize=12, verticalalignment='top')  # Text for number of monomers
ax2.set_xlim(0, time_steps)  # Set x limits for the plot
ax2.set_ylim(-0.5, 40)  # Set y limits for the plot
ax2.set_title("Number of Monomers and Theoretical Fit")  # Title for the monomer plot
ax2.set_xlabel("Time")  # Label for x-axis
ax2.set_ylabel("Number of Monomers")  # Label for y-axis
ax2.legend()  # Show legend for the plot

# Plot for the number of islands
line2, = ax3.plot([], [], color='orangered', label='Simulated Islands', alpha=0.9)  # Line for simulated islands
ax3.plot(1, -100, '--k', label='Theoretical Fit', lw=2)  # Placeholder for theoretical fit line
island_data, = ax3.plot(independent_island_counts, color='red', lw=2)  # Line for island data
island_text = ax3.text(0.6, 0.77, '', transform=ax3.transAxes, fontsize=12, verticalalignment='top')  # Text for number of islands
percent_text = ax3.text(0.42, 0.7, '', transform=ax3.transAxes, fontsize=12, verticalalignment='top')  # Text for island area percentage
ax3.set_xlim(0, time_steps)  # Set x limits for the plot
ax3.set_ylim(-0.5, 140)  # Set y limits for the plot
ax3.set_title("Number of Islands and Occupied Area")  # Title for the island plot
ax3.set_xlabel("Time")  # Label for x-axis
ax3.set_ylabel("Number of Islands")  # Label for y-axis
ax3.legend()  # Show legend for the plot

# Add theoretical curves with black dashed lines for n1 and n
fit_line_n1, = ax2.plot([], [], '--k', lw=2)  # Fitting line for n1
fit_line_n, = ax3.plot([], [], '--k', lw=2)  # Fitting line for n

# Text displays for the fitted parameters
d1_text = ax2.text(0.7, 0.67, '', transform=ax2.transAxes, fontsize=12, verticalalignment='top')  # Text for D1
a_text = ax2.text(0.7, 0.57, '', transform=ax2.transAxes, fontsize=12, verticalalignment='top')  # Text for a

# Function to update the animation
def animate(frame):
    update(frame)  # Call the update function for the current frame
    
    # Update the plots for the number of monomers and independent islands
    line1.set_data(np.arange(len(monomer_counts)), monomer_counts)  # Update monomer line data
    line2.set_data(np.arange(len(independent_island_counts)), independent_island_counts)  # Update island line data
    
    # Update text with current counts
    monomer_text.set_text(f"Monomers: {monomer_counts[-1]}")  # Display the current number of monomers
    island_text.set_text(f"Independent Islands: {independent_island_counts[-1]}")  # Display the current number of islands
    percent_text.set_text(f"Occupied Area by Islands: {island_percentages[-1]:.2f}%")  # Display the percentage of area occupied by islands
    
    # Optionally: adjust y-axis limits for better visualization
    ax2.set_ylim(0, max(max(monomer_counts), 40))  # Adjust upper limit for the number of monomers
    ax3.set_ylim(0, max(max(independent_island_counts), 140))  # Adjust upper limit for the number of islands
    
    return [mat, line1, line2, monomer_text, island_text, percent_text, d1_text, a_text]  # Return updated elements for animation

# Create the animation
ani = animation.FuncAnimation(fig, animate, frames=range(time_steps), repeat=False, interval=1)  # Animate the update function


#writer = animation.FFMpegWriter(fps=25)  # Create a video writer object
#ani.save('Pacman_Sim.mp4', writer=writer)  # Save the animation as an MP4 video

plt.tight_layout()  # Adjust layout to fit elements neatly
plt.show()  # Display the plots and animation

