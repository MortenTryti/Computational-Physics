import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import sys
#sns.set_style("darkgrid")


#Loading in the data from C++
complex_mat = pa.cx_cube()
prob_mat = pa.cube()

complex_mat.load("data/Task8_imag_cube.bin")

prob_mat = complex_mat @ pa.conj(complex_mat)


complex_mat = np.array(complex_mat)
prob_mat = np.array(prob_mat)
prob_mat = np.real(prob_mat)


real_mat = np.real(complex_mat)
imag_real = np.imag(complex_mat)

timelist = ["t=0","t=0.001","t=0.002"]

#Plotting the probability
for i in range(3):
    fig = plt.figure(i)
    im = plt.imshow(prob_mat[i], extent=[0,1,0,1])
    cbar = plt.colorbar(im)
    cbar.set_label(r"$|\psi|^2$ [-]")
    plt.xlabel("x [-]")
    plt.ylabel("y [-]")    
    plt.savefig("plots/task8_prob_"+str(i)+".pdf")

#Plotting the real
for i in range(3):
    plt.figure(i+3)
    im = plt.imshow(np.real(complex_mat[i]), extent=[0,1,0,1])
    cbar = plt.colorbar(im)
    cbar.set_label(r"Re($\psi$) [-]")
    plt.xlabel("x [-]")
    plt.ylabel("y [-]")    
    plt.savefig("plots/task8_real_"+str(i)+".pdf")

#Plotting the imaginary
for i in range(3):
    plt.figure(i+6)
    im = plt.imshow(np.imag(complex_mat[i]), extent=[0,1,0,1])
    cbar = plt.colorbar(im)
    cbar.set_label(r"Im($\psi$) [-]")
    plt.xlabel("x [-]")
    plt.ylabel("y [-]")    
    plt.savefig("plots/task8_imag_"+str(i)+".pdf")
plt.show()

