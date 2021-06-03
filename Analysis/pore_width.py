import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import hole2
import matplotlib.pyplot as plt
import os


variant= "M29"


#in pymol, select three atoms that line the pore and console: print(cmd.get_coords('sele', 1))
xyz = [[123.248,  73.546,  55.547],
 [110.482,  59.22,   56.434],
 [120.664,  55.43,   56.861]]

coordxyz = np.mean(xyz, axis = 0)
vec1 = np.subtract(xyz[0],xyz[1])
vec2 = np.subtract(xyz[0],xyz[2])
vectxyz = np.cross(vec1,vec2)
print(coordxyz, vectxyz)


u= mda.Universe( "../"+variant+".prmtop", "../"+variant+"_combined_imaged_sk10_r1.nc","../"+variant+"_combined_imaged_sk10_r2.nc","../"+variant+"_combined_imaged_sk10_r3.nc")
ha = hole2.HoleAnalysis(u, executable='/usr/local/hole2/exe/hole', end_radius=10.0, cpoint=coordxyz, cvect=vectxyz)
ha.run(1,1000,10)

plot=ha.plot_mean_profile(bins=100)
meanprof=plot.get_figure()
axls=meanprof.get_axes()
meanprof.savefig(variant+'profile.png')


minrad=ha.min_radius()
minrad.tofile('../'+variant+' minrad.txt', sep=",", format="%s")
x, y = minrad[:, 0], minrad[:, 1]

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
fs=12/2
axs[0].plot(x,y, color='black')
#axs[0].set_title('Minimum pore diameter over time', fontsize=fs)
axs[0].set_ylabel('HOLE radius (Å)', fontsize=fs-2)
violin_parts = axs[1].violinplot(y,showmeans=True, showextrema=True, showmedians=False)
for pc in violin_parts['bodies']:
    pc.set_facecolor('lightgrey')
violin_parts['cbars'].set_alpha(0)
violin_parts['cmeans'].set_color('black')
violin_parts['cmaxes'].set_color('black')
violin_parts['cmins'].set_color('black') 
#axs[1].set_title('Minimum pore diameter', fontsize=fs)
#axs[1].set_ylabel('HOLE radius (Å)', fontsize=fs-2)
axs[1].set_xticklabels([])
axs[1].set_xticks([])
fig.subplots_adjust(wspace=0.2)
fig.savefig(variant+'width.pdf')
plt.show()