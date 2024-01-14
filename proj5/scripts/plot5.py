import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# sns.set_style("darkgrid")


#Loading in the data from C++
S = [1,2,3]
for i in range(len(S)):
    u0_mat = pa.mat()
    u0_mat.load("data/M"+str(S[i])+".bin")
    u0_mat = np.array(u0_mat)
    plt.figure(i)
    im = plt.imshow(u0_mat, cmap='Greys', extent=[0,1,0,1])
    cbar = plt.colorbar(im)
    cbar.set_label(r"V [-]")
    plt.xlabel("x [-]")
    plt.ylabel("y [-]")    

    plt.savefig("plots/plot5_"+str(S[i])+".pdf")
plt.show()

