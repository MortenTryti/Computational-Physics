import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

method = ["RK4","Euler"]
for elem in method:

    #Importing all the positions and time, could be a 4 loop but Ctrl + C stronk
    tr1 = pa.mat()
    tr2 = pa.mat()
    tr3 = pa.mat()
    tr4 = pa.mat()

    tr1.load("TXYZcords"+elem+"false4000.bin")
    tr2.load("TXYZcords"+elem+"false8000.bin")
    tr3.load("TXYZcords"+elem+"false16000.bin")
    tr4.load("TXYZcords"+elem+"false32000.bin")
    tr1 = np.array(tr1)
    tr2 = np.array(tr2)
    tr3 = np.array(tr3)
    tr4 = np.array(tr4)

    t1,r1 = tr1[:,0], tr1[:,1:4]
    t2,r2 = tr2[:,0], tr2[:,1:4]
    t3,r3 = tr3[:,0], tr3[:,1:4]
    t4,r4 = tr4[:,0], tr4[:,1:4]

    T = [t1,t2,t3,t4]
    R = [r1,r2,r3,r4]

    #Stuff
    b0 = 9.648e1
    v0 = 2.41e6
    d = 500
    f = 0.1
    omega = 2
    q1 = 1
    m1 = 40.078
    wz2 = 2*q1*v0/(m1*d**2)
    w0 = q1*b0/m1
    wp = (w0+np.sqrt(w0**2-2*wz2))/2
    wm = (w0-np.sqrt(w0**2-2*wz2))/2
    Vel0 = 25
    Ap = (Vel0+wm*r1[0][0])/(wm-wp)
    Am =-(Vel0+wp*r1[0][0])/(wm-wp)

    #Exact solution to one particle problem
    #Returns the 3-vector r
    def r_exact(t,z0):
        f = Ap*np.exp(-1j*wp*t) + Am*np.exp(-1j*wm*t)
        r_x = np.real(f)
        r_y = np.imag(f)
        r_z = z0*np.cos(np.sqrt(wz2)*t)
        return np.array([r_x,r_y,r_z])


    #Define relative error function
    def rel_err(r_n,t):
        
        return np.linalg.norm(r_n-r_exact(t,r_n[0][2]).T,axis=1 )/np.linalg.norm(r_exact(t,r_n[0][2]),axis = 0)

    def delta(r_n,t):
        return np.max(np.linalg.norm(r_n-r_exact(t,r_n[0][2]).T,axis=1 ))

    # List of the number of datapoints
    counter = ["4000","8000","16000","32000"]

    for ri,ti,li in zip(R,T,counter):
        plt.semilogy(ti,rel_err(ri,ti),"-",label = "n="+li) 

    plt.legend()
    plt.xlabel("t [$\mu s$]")
    plt.ylabel("$log\epsilon_r$ [-]")
    plt.savefig(elem+"rel_err.pdf")
    plt.show()

    r_err = 0
    h_k = 50/np.array([4000,8000,16000,32000]) # mu s
    # Run over 2,3,4 which for lists would be 1,2,3
    for i in [1,2,3]:
        
        r_err += 1/3 * np.log(delta(R[i],T[i])/delta(R[i-1],T[i-1]))/np.log(h_k[i]/h_k[i-1]) 
    print("Rate of convergence is for "+elem )
    print(r_err)

