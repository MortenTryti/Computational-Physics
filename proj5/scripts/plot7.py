import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")


#Loading in the data from C++

textlist = ["noPot","Wpot"]
Label = ["now barrier","with double-slit barrier"]
for i in range(len(textlist)):
# for elem in textlist:
    u0_mat = pa.mat()


    u0_mat.load("data/probTask7"+textlist[i]+".bin")
    prob_vec = np.array(u0_mat)

    timevec = np.linspace(0,0.008,len(u0_mat))

    plt.ylabel(r"1-$|\psi|^2$ [-]")
    plt.xlabel("t [-]")
    plt.plot(timevec,prob_vec, label = Label[i])
    plt.yscale("log")
    plt.legend()
    plt.savefig("plots/task7"+textlist[i]+".pdf")
    plt.show()
