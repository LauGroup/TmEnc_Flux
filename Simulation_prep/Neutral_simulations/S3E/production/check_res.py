import mdtraj as md
import numpy as np 


def check_match(spec, match):
    resnames = [atom.residue.name for atom in traj.top.atoms if atom.residue.index in spec]
    resnames_set = set(resnames)
    print('set of resn')
    print(resnames_set)
    assert(len(resnames_set) == 1) #check only selected one residue
    assert(list(resnames_set)[0] == match) #check that it matches a specific res


traj = md.load("M26_npt.ncrst", top="M26.prmtop")
helix_start_spec = np.asarray([1,265,529,793,1057]) -1
helix_end_spec = np.asarray([121,385,649,913,1177]) -1
beta_start_spec = np.asarray([221,485,749,1013,1277]) -1
beta_end_spec = np.asarray([254,518,782,1046,1310]) -1


check_match(helix_start_spec, 'MET')
check_match(helix_end_spec, 'ARG')
check_match(beta_start_spec, 'ASP')
check_match(beta_end_spec, 'VAL')

print('done')
restraintmask = "'(:"
for i in range(5):
    restraintmask += f"{helix_start_spec[i]+1}-{helix_end_spec[i]+1},{beta_start_spec[i]+1}-{beta_end_spec[i]+1}"
    if i < 4:
        restraintmask += ","
restraintmask += ") & (@CA)'"
print(restraintmask)
