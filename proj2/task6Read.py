import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    #Datapoints
    N = [10,100]


    for i in N:
        J = np.loadtxt(f"task6Output{i}_Jacobi.txt", dtype = str)
        J = J.transpose()
        
        A = np.loadtxt(f"task6Output{i}_analytic.txt", dtype = str)
        A = A.transpose()
        

        x = A[1]
        print(len(x))


        a1, a2, a3 = A[3], A[5], A[7]
        j1, j2, j3 = J[3], J[5], J[7]
        # plt.figure(i)
        plt.figure(i,figsize=(10,6))
        plt.plot(x.astype(float),j1.astype(float), "--", label = f"N={i}, Jacobi, 1st")
        plt.plot(x.astype(float),j2.astype(float), "--", label = f"N={i}, Jacobi, 2nd")
        plt.plot(x.astype(float),j3.astype(float), "--", label = f"N={i}, Jacobi, 3rd")
        plt.plot(x.astype(float),a1.astype(float), "*", label = f"N={i}, analytic, 1st")
        plt.plot(x.astype(float),a2.astype(float), "*", label = f"N={i}, analytic, 2nd")
        plt.plot(x.astype(float),a3.astype(float), "*", label = f"N={i}, analytic, 3rd")        
        plt.legend()
        plt.xlabel("x [-]")
        plt.ylabel(" y [-]")
        plt.grid()
        plt.savefig("task6Plot_"+str(i)+".pdf")
        plt.show()

        

