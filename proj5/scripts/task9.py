import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")


#Loading in the data from C++
namelist = ["oneSlit","twoSlit","threeSlit"]
Label = ["one slit","two slits","three slits"]

# for name in namelist:
for i in range(len(namelist)):

    wavefunc_mat = pa.cx_mat()
    wavefunc_mat.load("data/Task9_"+namelist[i]+"_vec.bin")

    prob_mat = wavefunc_mat @ pa.conj(wavefunc_mat)


    prob_mat = np.array(prob_mat)
    prob_mat = np.real(prob_mat)

    x = np.linspace(0,1,len(prob_mat))
    plt.figure(i)
    plt.plot(x,prob_mat,label=Label[i])
    plt.xlabel("x [-]")
    plt.ylabel(r"$|\psi|^2$ [-]")
    plt.legend()
    plt.savefig("plots/task9_"+namelist[i]+".pdf")
plt.show()
