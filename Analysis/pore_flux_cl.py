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
cyl_rad = 7.5 #angstrom
cyl_height = 10 #angstrom
twiddle_z_pore = -3.0 #shift pore down slightly

print("ENCAPSULIN_FLUX_COUNTER\n\n")
mda.start_logging()
u = mda.Universe("no_water.M29_ions_1264.prmtop", "M29_ions_1264_combined_imaged_full_no_water.nc")

#select all the F180s
center_pore_ag = u.select_atoms("(resid 180 437 694 951 1208) and name CA")
assert(len(set(center_pore_ag.resnames)) == 1)
assert(list(set(center_pore_ag.resnames))[0] == "PHE")

#compute minimum distance between F180 Calphas  and set the cylinder radius to this value
self_distances = distances.self_distance_array(center_pore_ag.positions)
max_dist = np.max(self_distances)
print(f"Max distance between F180s is {max_dist}")

center_pore_com = center_pore_ag.center_of_mass()
print(f"Pore COM is  {[center_pore_com[0], center_pore_com[1], center_pore_com[2] + twiddle_z_pore]}\nSetting cyl_center to this value\n")

# select all the ions
ion_ag = u.select_atoms("resname Cl-")
# check that ions are atomic
assert(ion_ag.n_atoms == ion_ag.n_residues)
print(f"Number of ions is {len(ion_ag)}\n")
print(f"Number of frames is  {len(u.trajectory)}\n")

result = np.zeros([len(u.trajectory), ion_ag.n_atoms])

for i, ts in enumerate(ProgressBar(u.trajectory)):
    center_pore_com_active = np.asarray([center_pore_com[0], center_pore_com[1], center_pore_com[2] + twiddle_z_pore])
    result[i,:] = check_in_cylinder_upper_or_lower(ion_ag.positions, center_pore_com_active, cyl_rad, cyl_height)

print("finished main loop\n")

print("## RESULTS ##\n")
transitions = count_transitions(result)
print(f"transitions with oscillations {transitions}\n")

transitions, total_events = count_transitions_exit(result)
print(f"transitions discounting oscillations {transitions}\n")
print(f"total pore interaction events {total_events}\n")

plt.matshow(result, cmap=matplotlib.cm.bwr, interpolation="none", aspect="auto")
plt.show()


    




