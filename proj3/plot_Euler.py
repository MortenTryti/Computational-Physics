import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
# pa.load("zcord.bin")

C = pa.mat()
C.load("Zcord_Euler.bin")
C = np.array(C)
print(C[:,0])

plt.plot(C[:,1], C[:,0])
plt.savefig("Zcord_interaction_Euler.pdf")
plt.show()
