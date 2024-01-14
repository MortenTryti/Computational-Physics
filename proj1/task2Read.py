import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("task2Output.txt", dtype = str)

x = x.transpose()


xn = x[1]
un = x[3]


for elem1,elem2 in zip(xn,un):
    plt.plot(float(elem1),float(elem2), "+r")
    


plt.grid()
plt.xlabel("x [-]")
plt.ylabel("u(x) [-]")
plt.savefig("task2.pdf")
plt.show()