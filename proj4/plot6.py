import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

#Loading in the data from C++
eps_mat = pa.mat()
eps_mat.load("eps.bin")
eps_mat = np.array(eps_mat)


# Splitting into what i want
epsT1 = eps_mat[:,0]#[0::4]
epsT2 = eps_mat[:,1]#[0::4]

#Checks to see which elements are in the lists, was used to analyze app. binwidth etc
list1 = []
list2 = []
for elem in epsT1:
    if elem not in list1:
        list1.append(elem)

for elem in epsT2:
    if elem not in list2:
        list2.append(elem)



sortlist1 = []
sortlist2 = []
for elem in np.sort(np.array(list1)):
    sortlist1.append(elem)
for elem in np.sort(np.array(list2)):
    sortlist2.append(elem)

#A heavy way of doing this, should use min and max
bwidth1 = np.arange(sortlist1[0],sortlist1[-1],0.01)
bwidth2 = np.arange(sortlist2[0],sortlist2[-1],0.01)


weights1 = 1/len(epsT1)
weights2 = 1/len(epsT2) 

weightlist1 = [weights1 for i in range(len(epsT1))]
weightlist2 = [weights2 for i in range(len(epsT2))]


#This serves no purpose now, but I dont want to delete it 
counts1, bins1 = np.histogram(epsT1,bins = bwidth1,weights=weightlist1)
counts2, bins2 = np.histogram(epsT2,bins = bwidth2,weights=weightlist2)

# Checking to see if the sum of the bins add up to one, normalisation bby
#Density option would only give percentage, not decimals
#print(np.sum(counts1),np.sum(counts2))



print(f"The variance of eps at T=1 is var(eps) = {np.var(epsT1)}")
print(f"The variance of eps at T=2.4 is var(eps) = {np.var(epsT2)}")

#Plotting!
plt.hist(epsT1,bins = bwidth1, weights=weightlist1)
plt.xlabel("$\epsilon$ [J]")
plt.ylabel("$p_\epsilon$ [-]")
plt.legend(["T=1"])
plt.savefig("Task6T1.pdf")
plt.show()
plt.hist(epsT2,bwidth2,weights=weightlist2)
plt.xlabel("$\epsilon$ [J]")
plt.ylabel("$p_\epsilon$ [-]")
plt.legend(["T=2.4"])
plt.savefig("Task6T24.pdf")
plt.show()



