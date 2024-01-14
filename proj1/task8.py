import numpy as np
import matplotlib.pyplot as plt
from task7Read import u
from tabulate import tabulate


def relErr(u,v):
    delta = np.absolute((u-v)/u)

    return np.log10(delta)


def totErr(u,v):
    delta = np.absolute(u-v)
    return np.log10(delta)


N = [10,100,1000,10000]

#a
for i in N:
    #loading data and formating it to plotting.
    a = np.loadtxt(f"task7Output{i}.txt", dtype = str)
    a = a.transpose()
    

    x,v = a[1].astype(float), a[3].astype(float)
    #cutting off endpoints
    x,v = x[1:-1], v[1:-1]
    uarr = u(x)
    #Defining total error
    delTot = totErr(uarr,v)
    #plotting
    plt.plot(x,delTot, label =f"N={i}")


plt.xlabel("x  [-]")
plt.ylabel("$\log_{10}(u-v)$   [-]")
plt.title("Total error between numerical and analytic solution")
plt.grid()
plt.legend()
plt.savefig("task8_totErr.pdf")
plt.show()

#b
for i in N:
    #loading data and formating it to plotting
    a = np.loadtxt(f"task7Output{i}.txt", dtype = str)
    a = a.transpose()
    

    x,v = a[1].astype(float), a[3].astype(float)
    x,v = x[1:-1], v[1:-1]
    uarr = u(x)
    
    #Defining relative error
    delRel = relErr(uarr,v)
    plt.plot(x,delRel, label =f"N={i}")

plt.xlabel("x  [-]")
plt.ylabel(r"$\log_{10}(\frac{u-v}{u})$   [-]")
plt.title("Relative error between numerical  and analytic solution")
plt.grid()
plt.legend()
plt.savefig("task8_relErr.pdf")
plt.show()

#c 
N = [10,100,1000,10000,100000,1000000,10000000]
data = []

for i in N:
    #loading data and formating it to plotting
    a = np.loadtxt(f"task7Output{i}.txt", dtype = str)
    a = a.transpose()

    x,v = a[1].astype(float), a[3].astype(float)
    x,v = x[1:-1], v[1:-1]
    ufunc = u(x)
    mmax = np.max(np.absolute((ufunc-v)/ufunc))
    data.append([i,np.log10(mmax)])
print(tabulate(data,headers =["Number of iterations", "Maximum relative error"]))
