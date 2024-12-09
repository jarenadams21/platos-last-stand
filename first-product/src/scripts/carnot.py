import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("carnot_data.csv")

# P-V diagram
plt.figure()
plt.plot(df["V"], df["P"], 'o-', label="Carnot cycle")
plt.xlabel("Volume")
plt.ylabel("Pressure")
plt.title("P-V Diagram of Carnot Cycle")
plt.legend()
plt.show()
