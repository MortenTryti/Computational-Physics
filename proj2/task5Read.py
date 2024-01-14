import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# fit function for estematic the transformation scale
def func(x, a, b, c):
    return a*x**b + c


if __name__ == "__main__":
    
    #Datapoints for fit 
    x2 = np.linspace(0,100,1000)

    a = np.loadtxt(f"task5Output.txt", dtype = str)
    a = a.transpose()
    x,y = a[1], a[3]

    x1 = np.log10(float(x[-1]))
    x0 = np.log10(float(x[0]))


    y1 = np.log10(float(y[-1]))
    y0 = np.log10(float(y[0]))
    
    a = (y1-y0)/(x1-x0)
    print("slope ="+str(a))

    # power fit of f(x) = ax^b + c
    initial_guess = [1] * 3  # Initial guess for all parameters
    fit_params, covariance = curve_fit(func, x, y, p0=initial_guess)

    plt.figure(1)
    plt.plot(x.astype(float),y.astype(float), "--", label = f"Number of transformations")
    plt.plot(x2, func(x2, *fit_params), color='red', linestyle='solid', label=f"fit f(N) = {fit_params[0].round(3)}N^{fit_params[1].round(3)}  {fit_params[2].round(3)}")
    plt.legend()
    plt.xlabel("Matrix size N")
    plt.ylabel("Number of transformations")
    plt.grid()
    plt.savefig("task5Plot_fit.pdf")

    plt.figure(2)
    plt.axis()
    plt.plot(x.astype(float),y.astype(float), "--", label = f"Number of transformations")
    plt.xscale("log", base=10)
    plt.yscale("log", base=10)    
    plt.legend()
    plt.xlabel("Matrix size N")
    plt.ylabel("Number of transformations")
    plt.grid()
    plt.grid(which="major",alpha=0.6)
    plt.grid(which="minor",alpha=0.6)    
    plt.savefig("task5Plot_log.pdf")
    
    plt.show()

    print(f"Fit parameters: {fit_params}")

