import numpy as np
import matplotlib.pyplot as plt

# Analytic function
def u(x):
    return 1 - (1-np.exp(-10))*x - np.exp(-10*x)

if __name__ == "__main__":
    #Datapoints
    x =np.linspace(0,1,1000)
    plt.plot(x,u(x), "-b", label = "u(x)" )

    N = [10,100,1000,10000]


    for i in N:
        a = np.loadtxt(f"task7Output{i}.txt", dtype = str)
        a = a.transpose()
        

        x,v = a[1], a[3]

        plt.plot(x.astype(float),v.astype(float), "--", label = f"N={i}")
    plt.legend()
    plt.xlabel("x [-]")
    plt.ylabel(" u(x) [-]")
    plt.grid()
    plt.savefig("task7Plot.pdf")
    plt.show()

        