import mdtraj as md
import numpy as np


def check_match(spec, match):
    resnames = [
        atom.residue.name for atom in traj.top.atoms if atom.residue.index in spec]
    resnames_set = set(resnames)
    print('set of resn')
    print(resnames_set)
    assert(len(resnames_set) == 1)  # check only selected one residue
    # check that it matches a specific res
    assert(list(resnames_set)[0] == match)


traj = md.load("M8_npt.ncrst", top="M8.prmtop")
helix_start_spec = np.asarray([
1,
262,
523,
784,
1045,
                               
                               
                               
                               ]) - 1
helix_end_spec = np.asarray([121,
382,
643,
904,
1165,
                             
                             
                             
                             ]) - 1
beta_start_spec = np.asarray([217,
478,
739,
1000,
1261,
                              
                              
                              
                              ]) - 1
beta_end_spec = np.asarray([250,
511,
772,
1033,
1294,
                        
                            ]) - 1


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
restraintmask += ") & @CA'"
print(restraintmask)
