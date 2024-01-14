import sys
import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

# Load data
f = sys.argv[1]
C = pa.mat()
C.load("hundred_particles_"+str(f)+".bin")
C = np.array(C)
omega = C[:,0]
particles = C[:,1]
plt.plot(omega, particles, label="$f="+str(f)+"$")
plt.xlabel(r'$\omega_V$ $[MHz]$')
plt.ylabel('Fraction of particles still trapped $[-]$')
plt.legend()
plt.savefig("particles_inside_trap_"+str(f)+".pdf")
plt.show()
