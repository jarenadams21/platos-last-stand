import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the CSV data
# For demonstration, assume you've saved the posted CSV data into "qgp_data.csv"
df = pd.read_csv("qgp_data.csv")

# Parse the 'Era' column to distinguish radiation vs matter
# Let's assign a color based on era:
colors = df['Era'].map(lambda e: 'blue' if e.strip() == 'radiation' else 'red')

# Extract columns of interest
phi = df['phi'].astype(float)
pot = df['Potential'].astype(float)
kin = df['Kinetic'].astype(float)

# Because values vary wildly, consider using logarithms for better visual spread
# Avoid log of zero or negative by adding small offsets if needed:
import numpy as np

phi_log = np.log10(np.abs(phi) + 1e-60)  # add offset to avoid log(0)
pot_log = np.log10(np.abs(pot) + 1e-60)
kin_log = np.log10(np.abs(kin) + 1e-60)

# Create a 3D plot
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')

# Plot the data: each point (phi_log, pot_log, kin_log) with color by era
sc = ax.scatter(phi_log, pot_log, kin_log, c=colors, s=5, alpha=0.7)

# Add labels
ax.set_xlabel('log10(|phi|)')
ax.set_ylabel('log10(|Potential|)')
ax.set_zlabel('log10(|Kinetic|)')
ax.set_title('Quantum Gravity Plasma Field Evolution')

# Add a legend: radiation (blue), matter (red)
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Radiation Era', markerfacecolor='blue', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Matter Era', markerfacecolor='red', markersize=8)
]
ax.legend(handles=legend_elements, loc='best')

plt.tight_layout()
plt.show()
