import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")


#Loading in the data from C++
u0_mat = pa.mat()
u0_mat.load("../data/u0_prob_sq.bin")
u0_mat = np.array(u0_mat)



sns.heatmap(u0_mat)
plt.show()
