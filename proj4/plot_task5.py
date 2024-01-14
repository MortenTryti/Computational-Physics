import sys
import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

# Load data
# f = sys.argv[1]

vertical = np.linspace(0,2e4,5,dtype=int)

dict = {
  True: "Ordered",
  False: "Unordered"}
for elem in [True,False]:
    C = pa.mat()
    if (elem == True):
        C.load("task5True.bin")
    if (elem ==False):
        C.load("task5false.bin")
    C = np.array(C)

    cycles = C[:,0]#[0::100]
    eps1_cumsum = C[:,1]#[0::100]
    m1_cumsum = C[:,2]#[0::100]
    eps2_cumsum = C[:,3]#[0::100]
    m2_cumsum = C[:,4]#[0::100]


    plt.figure(1)
    
    plt.xlabel('cycles [-]')
    plt.ylabel(r'<e>[ J ]')
    if (elem==True):
        for elem2 in vertical:
            plt.axvline(elem2,linestyle = "--",label = f"cycles={elem2}")
    plt.plot(cycles, eps1_cumsum, label="T=1: Config = "+dict[elem])
    plt.legend()
    plt.savefig("Task5epsT1.pdf")
    plt.figure(2)
    if (elem==True):
        for elem2 in vertical:
            plt.axvline(elem2,linestyle = "--",label = f"cycles={elem2}")
    plt.plot(cycles, m1_cumsum, label="T=1: Config = "+dict[elem])
    plt.xlabel('cycles [-]')
    plt.ylabel(r'<|m|>[-]')
    
    plt.legend()
    plt.savefig("Task5mT1.pdf")
    plt.figure(3)
    if (elem==True):
        for elem2 in vertical:
            plt.axvline(elem2,linestyle = "--",label = f"cycles={elem2}")
    plt.plot(cycles, eps2_cumsum, label="T=2.4: Config = "+dict[elem])
    plt.xlabel('cycles [-]')
    plt.ylabel(r'<e>[ J ]')
    plt.legend()
    
    plt.savefig("Task5epsT2.pdf")

    plt.figure(4)
    if (elem==True):
        for elem2 in vertical:
            plt.axvline(elem2,linestyle = "--",label = f"cycles={elem2}")
    plt.plot(cycles, m2_cumsum, label="T=2.4: Config = "+dict[elem])
    plt.xlabel('cycles [-]')
    plt.ylabel(r'<|m|>[-]')
    
    plt.legend()
    plt.savefig("Task5mT2.pdf")

plt.show()
