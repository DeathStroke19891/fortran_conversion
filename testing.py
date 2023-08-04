import numpy as np

def kubo_conductivity(eigenvalues, eigenvectors, band_hopping, chemical_potential, temperature, energy_window):
    num_bands, num_sites = eigenvalues.shape

    # Calculate density of states
    dos = np.zeros_like(eigenvalues)
    for i in range(num_bands):
        for j in range(num_sites):
            dos[i, j] = np.sum(np.exp(-(eigenvalues[i, j] - chemical_potential) / (temperature * 8.617333262145e-5))) / temperature

    # Calculate site occupancy (0 for singly occupied, 1 for doubly occupied)
    site_occupancy = np.where(dos > 1, 1, 0)

    # Calculate dc conductivity
    conductivity = 0
    for i in range(num_bands):
        for j in range(num_sites):
            for k in range(num_bands):
                for l in range(num_sites):
                    if abs(eigenvalues[i, j] - eigenvalues[k, l]) < energy_window:
                        conductivity += (1 - site_occupancy[i, j]) * site_occupancy[k, l] * \
                                       np.abs(np.dot(eigenvectors[i, j], np.dot(band_hopping, eigenvectors[k, l])))**2

    return conductivity

# Example input parameters
num_bands = 3
num_sites = 4
eigenvalues = np.random.rand(num_bands, num_sites)
eigenvectors = np.random.rand(num_bands, num_sites, num_sites)
band_hopping = np.random.rand(num_sites, num_sites)
chemical_potential = 0.5
temperature = 300
energy_window = 0.1
print(eigenvalues)
print(eigenvectors)
print(band_hopping)

conductivity = kubo_conductivity(eigenvalues, eigenvectors, band_hopping, chemical_potential, temperature, energy_window)
print("DC Conductivity:", conductivity)
