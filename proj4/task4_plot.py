import sys
import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

C = pa.mat()
C.load("MC_Cycles.bin")
C = np.array(C)

print(C)
Tol = C[0,:].astype("float64")
typelist = ["<$\epsilon$>","<|m|>","$C_V$","$X$"]

for i in range(len(C[1::])):
    plt.plot(C[i+1,:],Tol,"o--",label = typelist[i])

plt.ylabel("tol")
plt.xlabel("# MC cycles")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.savefig("Task4C.pdf")
plt.show()


