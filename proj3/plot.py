import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")
# pa.load("zcord.bin")

C = pa.mat()
C.load("ZcordRK4.bin")
C = np.array(C)


plt.plot(C[:,1], C[:,0])
plt.xlabel("t [$\mu s$]")
plt.ylabel("z [$\mu m$]")
plt.show()
