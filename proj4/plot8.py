import ast
import sys
import pyarma as pa
import numpy as np
import pandas as pds
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import seaborn as sns
sns.set_style("darkgrid")



zoom = sys.argv[1]
T_points = sys.argv[2]
L = ast.literal_eval(sys.argv[3])

variables = pa.mat()
T_vec = pa.mat()

if int(zoom) == 0:
    T_vec.load("T_unzoom_"+str(T_points)+"_vec.bin")
elif int(zoom) == 1:
    T_vec.load("T_zoom_"+str(T_points)+"_vec.bin")    

T = np.array(T_vec)[:,0]


# number of varaibles to plot
n_var = 4
var_labels = [r'$ \langle \varepsilon \rangle $ [J]', r'$\langle |m| \rangle $ [-]', r'$C_V/k_B $ [-]', r'$\chi $ [1/J]']
var_names = ["epsilon", "mag", "CV", "chi"]
# colors
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

# vector to store the critical temperatures
Critical_T_CV = np.zeros(len(L))
Critical_T_chi = np.zeros(len(L))


for i in range(n_var):
    plt.figure(i)    
    for j in range(len(L)):
        if int(zoom) == 0:
            variables.load("variables_unzoom_"+str(T_points)+"_"+str(L[j])+".bin")            
        elif int(zoom) == 1:
            variables.load("variables_zoom_"+str(T_points)+"_"+str(L[j])+".bin")

        var = variables[:,i]
        plt.plot(T, var, 'o',label = "L = "+str(L[j]), color=colors[j])
        f = interpolate.interp1d(T, var, kind = "quadratic")
        T_new = np.linspace(T[0], T[-1], 1000)
        # find the critical temperature as the maximal of the interpolated function
        if int(zoom) == 1:
            if (var_names[i] == "CV"):
                max_element = pds.Series(f(T_new)).idxmax()
                Critical_T_CV[j] = T_new[max_element]
            elif (var_names[i] == "chi"):
                max_element = pds.Series(f(T_new)).idxmax()
                Critical_T_chi[j] = T_new[max_element]
        plt.plot(T_new, f(T_new), color=colors[j])
        plt.xlabel(r'T [J/$k_B$]')
        plt.ylabel(var_labels[i])
        plt.legend()
        if int(zoom) == 0:
            plt.savefig("Task8_unzoom_"+str(T_points)+"_"+var_names[i]+".pdf")
        elif int(zoom) == 1:
            plt.savefig("Task8_zoom_"+str(T_points)+"_"+var_names[i]+".pdf")

print(Critical_T_CV)
print(Critical_T_chi)

# make linear regression of the critical temperatures
if int(zoom) == 1:
    L_invers = [1/int(L[i]) for i in range(len(L))]
    L_invers_tick = ["1/"+str(L[i]) for i in range(len(L))]
    L_invers_tick_2 = [r'$\frac{1}{'+str(L[i])+'}$' for i in range(len(L))]    
    print(L_invers)
    # make a linear regression from the critical temperatures
    slope, intercept, r, p, se = stats.linregress(L_invers, Critical_T_CV)
    res_CV = stats.linregress(L_invers, Critical_T_CV)
    res_chi = stats.linregress(L_invers, Critical_T_chi)        

    a_CV = res_CV.slope
    a_chi = res_chi.slope
    a_chi_err = res_chi.stderr
    a_CV_err = res_CV.stderr
    print("a_CV = ", a_CV)
    print("a_chi = ", a_chi)
    print("a_CV_err = ", a_CV_err)
    print("a_chi_err = ", a_chi_err)

    Tc_inf_CV = res_CV.intercept
    Tc_inf_chi = res_chi.intercept
    Tc_inf_CV_err = res_CV.intercept_stderr
    Tc_inf_chi_err = res_chi.intercept_stderr
    print("Tc_inf_CV = ", Tc_inf_CV)
    print("Tc_inf_chi = ", Tc_inf_chi)
    print("Tc_inf_CV_err = ", Tc_inf_CV_err)
    print("Tc_inf_chi_err = ", Tc_inf_chi_err)
    
    L_fit = np.linspace(-1/1000, 1/30, 100)
    T_fit_CV = [a_CV*L_fit[i] + Tc_inf_CV for i in range(len(L_fit))]
    T_fit_chi = [a_chi*L_fit[i] + Tc_inf_chi for i in range(len(L_fit))]    
    
    print(res_CV.intercept)
    print(res_chi.intercept)    

    plt.figure(6)
    plt.ylim(2.25, 2.35)

    
    plt.plot(L_invers, Critical_T_CV, 'o',color=colors[0], label=r'$C_V$ data, $T_{C}(L)$')
    fit_lable_CV = r'$C_V$ fit,  $T_{C}(L) \approx $'+str(a_CV.round(3))+'$L^{-1} +$ '+str(Tc_inf_CV.round(3))
    plt.plot(L_fit, T_fit_CV,color=colors[0], label=fit_lable_CV)
    
    plt.plot(L_invers, Critical_T_chi, 'o',color=colors[1], label=r'$\chi$ data, $T_{C}(L)$')
    fit_lable_chi = r'$\chi$ fit, $T_{C}(L) \approx $'+str(a_chi.round(3))+'$L^{-1} +$ '+str(Tc_inf_chi.round(3))
    plt.plot(L_fit, T_fit_chi,color=colors[1], label=fit_lable_chi)

    plt.xlabel(r'$L^{-1}$ [-]')
    plt.ylabel(r'$T_c(L)$ [J/$k_B$]')
    # make new ticks on the x-axis as 1/L
    L_invers.insert(0,0)
    L_invers_tick_2.insert(0,"0")    
    plt.xticks(L_invers, L_invers_tick_2)
    plt.axhline(y=2.269, color='k',linestyle='--',label=r'theoretical, $T_C(L=\infty)\approx 2.269$')
    plt.legend()
    plt.savefig("Task8_regression_"+str(T_points)+"_critical_T.pdf")
plt.show()


