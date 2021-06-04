import MDAnalysis as mda
import MDAnalysis.analysis.rms
import matplotlib.pyplot as plt

variant = "M29"


u= mda.Universe( "../"+variant+".prmtop", "../"+variant+"_combined_imaged_sk10_r1.nc","../"+variant+"_combined_imaged_sk10_r2.nc","../"+variant+"_combined_imaged_sk10_r3.nc")
ref= mda.Universe( "../"+variant+".prmtop", "../"+variant+"_combined_imaged_sk10_r1.nc")
R = mda.analysis.rms.RMSD(u, ref,
           select="backbone",             # superimpose on whole backbone of the whole protein
           groupselections=["backbone and (resid 181-189 or resid 438-446 or resid 695-703 or resid 952-960 or resid 1209-1217)",
                           "backbone",
                           "backbone and (resid 330-342 or resid 71-86)"]) 
R.run()




rmsd = R.rmsd.T 
print(rmsd)
time = rmsd[1]
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.plot(time, rmsd[4], 'k-',  label="all")
ax.plot(time, rmsd[3], 'r-', label="pore")
ax.plot(time, rmsd[5], 'b-', label="edge")
ax.legend(loc="best")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"RMSD ($\AA$)")
fig.savefig(variant+"rmsd.pdf")