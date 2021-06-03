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


traj = md.load("M2_npt.ncrst", top="M2.prmtop")
helix_start_spec = np.asarray([1,
                               259,
                               517,
                               775,
                               1033]) - 1
helix_end_spec = np.asarray([121,
                             379,
                             637,
                             895,
                             1153]) - 1
beta_start_spec = np.asarray([214,
                              472,
                              730,
                              988,
                              1246]) - 1
beta_end_spec = np.asarray([247,
                            505,
                            763,
                            1021,
                            1279]) - 1


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
