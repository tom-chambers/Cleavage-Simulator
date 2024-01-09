import numpy as np
import matplotlib.pyplot as plt

# Set initiation (kini), translocation (kf), and dissociation (koff) rates for each enzyme

enzyme1_kini = 1
enzyme1_kf = 0.8
enzyme1_koff = 0.00

enzyme2_kini = 1
enzyme2_kf = 0.8
enzyme2_koff = 0.00

# Set cleavage rate (kcut)

kcut = 1

# Set both subunits to bound but not translocating

enzyme1_translocating = False
enzyme2_translocating = False


# Function for Monte Carlo simulation

def monte_carlo_step(probability):
    return np.random.rand() < probability

# Function to simulate translocation - edit plasmid length and maximum number of steps per simulatuon here

def simulate_translocation(plasmid_length=2500, max_steps=100000):
    global enzyme1_translocating, enzyme2_translocating
    global enzyme1_position, enzyme2_position

    enzyme1_position = 1
    enzyme2_position = plasmid_length
    steps = 0
    collision_step = None
    cleavage_occurred = False

    while steps < max_steps and enzyme1_position != enzyme2_position:
        
        # Randomly select which enzyme to consider for initiation/translocation
        
        active_enzyme = np.random.choice([1, 2])

        if active_enzyme == 1:
            
            # Enzyme 1 initiation
            if not enzyme1_translocating and monte_carlo_step(enzyme1_kini):
                enzyme1_translocating = True

            # Enzyme 1 translocation
            if enzyme1_translocating:
                if monte_carlo_step(enzyme1_koff):
                    enzyme1_translocating = False
                    enzyme1_position = 1
                elif monte_carlo_step(enzyme1_kf):
                    enzyme1_position = (enzyme1_position) + 1

        elif active_enzyme == 2:
            
            # Enzyme 2 initiation
            if not enzyme2_translocating and monte_carlo_step(enzyme2_kini):
                enzyme2_translocating = True

            # Enzyme 2 translocation
            if enzyme2_translocating:
                if monte_carlo_step(enzyme2_koff):
                    enzyme2_translocating = False
                    enzyme2_position = plasmid_length
                elif monte_carlo_step(enzyme2_kf):
                    enzyme2_position = (enzyme2_position) -1

        # Print positions - produces a readout of all steps during the simulations
        #print(f"Step {steps + 1}: Enzyme 1 at position {enzyme1_position}, Enzyme 2 at position {enzyme2_position}")

        # Checking for subunit collision
        
        if enzyme1_position == enzyme2_position:
            collision_step = steps + 1
            
            # Print collision sites - produces a readout of all collisions during the simulation
            # print(f"Enzymes collide at position {enzyme1_position} after {collision_step} steps.")
            
            # Check for cleavage
            
            if np.random.rand() < kcut:
                cleavage_occurred = True
            
            # Print cleavage sites - produces a readout of all cleavage sites during the simulation
                # print(f"Cleavage occurred at position {enzyme1_position}. Simulation ends.")
                
                break
                
            else:
                
                # Print failed cleavage attempts - produces a readout of all failed attempts at cleavage during the simulation
                # print("Cleavage failed. Enzymes return to initial positions and attempt initiation.")
                enzyme1_translocating = False
                enzyme2_translocating = False
                enzyme1_position = 1
                enzyme2_position = plasmid_length

        steps += 1

    if not cleavage_occurred and steps == max_steps:
        
        # Print complete simulations without enzyme cleavage
        # print("Simulation reached maximum steps without cleavage.")
        collision_step = None

    return enzyme1_position, enzyme2_position, collision_step

# Run multiple simulations - record cleavage position and number of steps required for cleavage for each simulation

num_simulations = 100
all_cleavage_positions = []
collision_steps = []

for i in range(num_simulations):
  
    result = simulate_translocation()
    all_cleavage_positions.append(result[0])  
    collision_steps.append(result[2])

# Print results

print("\nResults:")
for i in range(num_simulations):
    if collision_steps[i] is not None:
        print(f"Simulation {i + 1} - cleavage position: {all_cleavage_positions[i]}, number of steps: {collision_steps[i]}")
    else:
        print(f"Simulation {i + 1} - no cleavage occurred.")

# Plot results 

plt.figure(figsize=(10, 6))
for i in range(num_simulations):
    if collision_steps[i] is not None:
        plt.plot(all_cleavage_positions[i], 0, 'gx')

plt.legend(['Cleavage'])
plt.xlabel('Position on plasmid')
plt.title('Simulation of TI restriction endonuclease translocation and cleavage')
plt.xlim(0, 2500)
plt.yticks([])

plt.show()