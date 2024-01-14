import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")


for elem in ["true","false"]:
    TXYZ = pa.mat()
    TXYZ.load("TXYZcordsRK4"+elem+"4000.bin")
    TXYZ = np.array(TXYZ)

    VXYZ = pa.mat()
    VXYZ.load("XYZ_VelRK4"+elem+"4000.bin")
    VXYZ = np.array(VXYZ)

    T = TXYZ[:,0]

    X1 = TXYZ[:,1]
    Y1 = TXYZ[:,2]
    Z1 = TXYZ[:,3]

    X2 = TXYZ[:,4]
    Y2 = TXYZ[:,5]
    Z2 = TXYZ[:,6]

    VX1 = VXYZ[:,0]
    VY1 = VXYZ[:,1]
    VZ1 = VXYZ[:,2]

    VX2 = VXYZ[:,3]
    VY2 = VXYZ[:,4]
    VZ2 = VXYZ[:,5]

    plt.figure(1)
    plt.scatter(T[0],Z1[0], color = "orange")
    #plt.scatter(T[0],Z2[0], color = "dodgerblue")
    plt.plot(T,Z1,color = "orange",label = "P1")
    #plt.plot(T,Z2,color = "dodgerblue",label = "P2")
    plt.xlabel("t [$\mu s$]")
    plt.ylabel("z [$\mu m$]")
    plt.axis('equal')
    plt.legend()
    plt.savefig("t_z_plot"+elem+".pdf")

    plt.figure(2)
    plt.scatter(X1[0],Y1[0],color = "orange")
    plt.scatter(X2[0],Y2[0],color = "dodgerblue")
    plt.plot(X1,Y1,color = "orange",label = "P1")
    plt.plot(X2,Y2,color = "dodgerblue",label = "P2")
    plt.xlabel("x [$\mu s$]")
    plt.ylabel("y [$\mu m$]")
    plt.axis('equal')
    plt.legend()
    plt.savefig("x_yplot"+elem+".pdf")

    plt.figure(3)
    plt.scatter(X1[0],VX1[0],color = "orange")
    plt.scatter(X2[0],VX2[0],color = "dodgerblue")
    plt.plot(X1,VX1,color = "orange",label = "P1")
    plt.plot(X2,VX2,color = "dodgerblue",label = "P2")
    plt.xlabel("x [$\mu m$]")
    plt.ylabel("$v_x$ [$\mu m/\mu s$]")
    plt.axis('equal')
    plt.legend()
    plt.savefig("x_vx_yplot"+elem+".pdf")

    plt.figure(4)
    plt.scatter(Z1[0],VZ1[0],color = "orange")
    plt.scatter(Z2[0],VZ2[0],color = "dodgerblue")
    plt.plot(Z1,VZ1,color = "orange",label = "P1")
    plt.plot(Z2,VZ2,color = "dodgerblue",label = "P2")
    plt.xlabel("z [$\mu m$]")
    plt.ylabel("$v_z$ [$\mu m/\mu s$]")
    plt.axis('equal')
    plt.legend()
    plt.savefig("z_vz_yplot"+elem+".pdf")

    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(X1[0],Y1[0],Z1[0],color = "orange")
    ax.scatter(X2[0],Y2[0],Z2[0],color = "dodgerblue")
    ax.plot(X1,Y1,Z1,color = "orange",label = "P1")
    ax.plot(X2,Y2,Z2,color = "dodgerblue",label = "P2")
    ax.set_xlabel("x [$\mu m$]")
    ax.set_ylabel("y [$\mu m$]")
    ax.set_zlabel("z [$\mu m$]")
    plt.legend()
    plt.savefig("3Dplot"+elem+".pdf")
    plt.show()

