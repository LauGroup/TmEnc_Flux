import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.lib.log import ProgressBar




## we are looking for transitions through the pore:
# 1. assign a cylinder whose center plane is the center of the pore
# 2. assign each ion a state
#       ->  1 if in top half of cylinder
#       -> -1 in bottom half of cylinder
#       ->  0 if not in cylinder at all
#
# all transitions of interest are then of the form 1 -> -1 or -1 -> 1.
# this covers the case where we travel through the PBC because must travel
# through 0 state intermediate , thereby being discounted


# non-vectorised so quite slow
def check_in_cylinder_upper_or_lower(pos_array, cylinder_center, cylinder_rad, cylinder_height):
    n_atoms = pos_array.shape[0]
    result = np.zeros(n_atoms)
    n_in_cyl = 0
    for i in range(n_atoms):
        # xdim
        if (pos_array[i,0] < cylinder_center[0] + cylinder_rad) and (pos_array[i,0] > cylinder_center[0] - cylinder_rad):
            # inside in x now check ydim 
            if (pos_array[i,1] < cylinder_center[1] + cylinder_rad) and (pos_array[i,1] > cylinder_center[1] - cylinder_rad):
                # inside in x and y now check if its in the upper or lower half of cylinder
                if (pos_array[i,2] < cylinder_center[2] + cylinder_height) and (pos_array[i,2] > cylinder_center[2]):
                    # in the upper half of the cylinder
                    result[i] = 1
                    n_in_cyl += 1
                elif (pos_array[i,2] > cylinder_center[2] - cylinder_height) and (pos_array[i,2] < cylinder_center[2]):
                    # in the lower half of the cylinder
                    result[i] = -1 
                    n_in_cyl += 1
                else:
                    # would set to 0 but is already
                    pass
    # print(n_in_cyl) 
    return result

#NOTE bad memory access pattern
def count_transitions(result_array):
    transitions = 0
    # slice out the timeseries for each residue
    for res in range(result_array.shape[1]):
        # make sure entry condition is right
        prev_val = result_array[0,res]
        for frame in range(result_array.shape[0]):
            val = result_array[frame,res]
            if val == prev_val:
                # no change, we arn't interested
                pass
            else:
                # we have either entered the cylinder or moved from one side to the other
                if (val == -1) and (prev_val == 1):
                    # this is a transition we are interested in -ve to +ve
                    transitions += 1
                elif (val == 1) and (prev_val == -1):
                    # this is a transition we are interested in +ve to -ve
                    transitions += 1
                else:
                    # either left or entered cylinder, we arn't intrested
                    pass
            prev_val = val

    return(transitions)


#NOTE bad memory access pattern
def count_transitions_exit(result_array):
    transitions = 0
    total_events = 0
    # slice out the timeseries for each residue
    for res in range(result_array.shape[1]):
        # make sure entry condition is right
        prev_val = result_array[0,res]
        prev_internal_flag = bool(result_array[0,res])
        if prev_internal_flag:
            cyl_entry =  result_array[0,res]       
        for frame in range(result_array.shape[0]):
            val = result_array[frame,res]
            internal_flag = bool(result_array[frame,res])
            
            # we have moved into the cylinder, track which side we came from
            if (not(prev_internal_flag) and internal_flag):
                cyl_entry = val
                total_events += 1
            
            if val == prev_val:
                # no change, we arn't interested
                pass

            else:
                # we were inside the cylinder on last iteration
                if prev_internal_flag:
                    if (val == 0) and (prev_val == 1):
                        # we have left the cylinder on the +ve side
                        if cyl_entry == -1:
                            # which is the opposite side to the one we came in 
                            transitions += 1
                    elif (val == 0) and (prev_val == -1):
                        # we have left the cylinder on the -ve side
                        if cyl_entry == 1:
                            # which is the opposite one to the one we came in on 
                            transitions += 1
                    else:
                        # we stayed inside the cylinder
                        pass
                else:
                    # we entered or left the cylinder without passing through
                    pass

            prev_internal_flag = bool(val)
            prev_val = val
        #print(f"res  {res}, trans {transitions}")

    return(transitions, total_events)



    


# backup constants 
cyl_rad = 15 #angstrom
cyl_height = 10 #angstrom
twiddle_z_pore = -0.0 #shift pore down slightly

print("ENCAPSULIN_FLUX_COUNTER\n\n")
mda.start_logging()
u = mda.Universe("no_water.M29_ions_1264.prmtop", "M29_ions_1264_combined_imaged_full_no_water.nc")

center_pore_ag_1 = u.select_atoms("(resid 25 97 506 99) and name CA")
print(center_pore_ag_1.resnames)
center_pore_ag_2 = u.select_atoms("(resid 796 868 870 249) and name CA")
print(center_pore_ag_2.resnames)
center_pore_ag_3 = u.select_atoms("(resid 1125 1127 1020 1053) and name CA")
print(center_pore_ag_3.resnames)
center_pore_ag_4 = u.select_atoms("(resid 611 613 1277 1193) and name CA")
print(center_pore_ag_4.resnames)
center_pore_ag_5 = u.select_atoms("(resid 354 356 763 282) and name CA")
print(center_pore_ag_5.resnames)


pore_ags = [center_pore_ag_1,center_pore_ag_2,center_pore_ag_3,center_pore_ag_4,center_pore_ag_5]#compute minimum distance between F180 Calphas  and set the cylinder radius to this value
for pore_ag in pore_ags:
    pore_ag_com = pore_ag.center_of_mass()
    print(f"Pore COM is  {[pore_ag_com[0], pore_ag_com[1], pore_ag_com[2] + twiddle_z_pore]}")

# select all the ions
ion_ag = u.select_atoms("resname Cl-")
# check that ions are atomic
assert(ion_ag.n_atoms == ion_ag.n_residues)
print(f"Number of ions is {len(ion_ag)}\n")
print(f"Number of frames is  {len(u.trajectory)}\n")

result = np.zeros([len(pore_ags),len(u.trajectory), ion_ag.n_atoms])

for i, ts in enumerate(ProgressBar(u.trajectory)):
    for j, pore_ag in enumerate(pore_ags):
        pore_ag_com = pore_ag.center_of_mass()
        center_pore_com_active = np.asarray([pore_ag_com[0], pore_ag_com[1], pore_ag_com[2] + twiddle_z_pore])
        result[j,i,:] = check_in_cylinder_upper_or_lower(ion_ag.positions, center_pore_com_active, cyl_rad, cyl_height)

print("finished main loop\n")
print("## RESULTS ##\n")
t_oscs = []
t_noosc = []
interactions = []

for i in range(len(pore_ags)):
    transitions = count_transitions(result[i,:,:])
    t_oscs.append(transitions)
    print(f"transitions with oscillations {transitions}\n")
    transitions, total_events = count_transitions_exit(result[i,:,:])
    t_noosc.append(transitions)
    interactions.append(total_events)
    print(f"transitions discounting oscillations {transitions}\n")
    print(f"total pore interaction events {total_events}\n")

arr = np.vstack([t_oscs,t_noosc, interactions])
print(arr.T)
np.savetxt("3xpores_cl.csv", arr.T, delimiter=",", fmt='%i')
# plt.matshow(result, cmap=matplotlib.cm.bwr, interpolation="none", aspect="auto")
# plt.show()


    




