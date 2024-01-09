CLEAVAGE SIMULATOR README

Cleavage Simulator enables the modelling of type I restriction endonuclease initiation, translocation, dissociation, and cleavage on a circular DNA plasmid. This software utilises a Monte Carlo simulation to model the stochastic nature of enzymes. Cleavage Simulator can run multiple simulations while retaining data relating to each individual event. 

This model assumes a single binding site for type I restriction endonucleases. However, the code could be modified to account for multiple binding sites. Additionally, the code could also be modified to model the mechanisms of other motor proteins that translocate nucleic acids, including helices and polymerases. It would also be possible to model the effect of drugs, such as small-molecule inhibitors, on motor proteins. 


Required libraries:

NumPy
Matplotlib


To run the software:

Edit the cleavage_simulator.py file to select desired parameters.

The initiation rate (kini), translocation rate (kf), and dissociation rate (koff) of each subunit can be adjusted (0 - 1):

enzyme1_kini = 0.4
enzyme1_kf = 1
enzyme1_koff = 0.002

enzyme2_kini = 0.4
enzyme2_kf = 1
enzyme2_koff = 0.002

The cleavage rate can be adjusted (0 - 1):

kcut = 1


The plasmid length and max steps can be adjusted within the simulate_translocation function:

def simulate_translocation(plasmid_length=2500, max_steps=100000)

The x-axis of the graphical output should be adjusted to match the plasmid length:

plt.xlim(0, 2500)


The total number of simulations can be adjusted:

num_simulations = 100

Run the software and a graphical output will be generated showing the locations of cleavage along with a readout of the result of each simulation (e.g., site of cleavage and how many steps, or no cleavage recorded).

By default, data relating to individual simulations (initiation, dissociation, number of steps, locations of collision, failed cleavage attempts) is not shown. This data can be shown by removing the # before the respective print functions:

# print(f"Step {steps + 1}: Enzyme 1 at position {enzyme1_position}, Enzyme 2 at position {enzyme2_position}")
# print(f"Enzymes collide at position {enzyme1_position} after {collision_step} steps.")
# print(f"Cleavage occurred at position {enzyme1_position}. Simulation ends.")
# print("Cleavage failed. Enzymes return to initial positions and attempt initiation.")
# print("Simulation reached maximum steps without cleavage.")
