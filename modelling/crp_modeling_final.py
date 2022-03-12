
##########################################################################################
############### IMP Modeling Script for Pinhead ###############
##########################################################################################


# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
import IMP.desmosome
import IMP.centrioles
import sys

runID = sys.argv[1]   # Specify the number of runs
run_output_dir = 'run_' + str(runID)


num_frames = 10000

max_temp = 4

# Identify data files
data_directory = "./data/"
gmm_data_D = data_directory+"gmm/EMD_9171_Pinhead_1-2_gmm.txt"                #Pinhead 1-2 only


# Topology File
topology_file = data_directory+"topology.txt"

##Paramaters
MPDBR_WEIGHT = 5             #Weight of MPDBR Restraint
Y2H_COIP_hetero_FACTOR = 2   #Factor to multiply the MPDBR weight to weigh Y2H/co-IP experimental evidence #Hetero protein interaction weighed < homo protein interaction
Y2H_COIP_homo_FACTOR = 4     #Factor to multiply the MPDBR weight to weigh Y2H/co-IP experimental evidence #Hetero protein interaction weighed < homo protein interaction
MPDBR_RES = 30               #Resolution of the representation to apply the MPDBR Restraint on
SPMGR_RES = 10               #Resolution of the representation to apply the SPMGR Restraint on
SPMGR_WEIGHT = 10            #Factor to multiply the SPMGR weight to weigh Y2H/co-IP experimental evidence for Mt binding
EM_WEIGHT = 10               #Factor to multiply the EM weight to weigh cryoEM experimental evidence for localization around MT regions

# FUNCTIONS AND WRAPPERS -----------------------------------------------
# wrapper for the MPDBR Restraint
class MinimumPairDistanceBindingRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist1, plist2, x0=0, kappa=1, label=None, weight=1):
        name = 'MinimumPairDistanceBindingRestraint%1%'
        super(MinimumPairDistanceBindingRestraint, self).__init__(model, name=name, label=label, weight=weight)
        l1 = IMP.container.ListSingletonContainer(mdl)
        l1.add(plist1)
        l2 = IMP.container.ListSingletonContainer(mdl)
        l2.add(plist2)
        bipartite_container = IMP.container.AllBipartitePairContainer(l1, l2)
        score = IMP.core.HarmonicSphereDistancePairScore(x0, kappa)
        res_main = IMP.container.MinimumPairRestraint(score, bipartite_container, 1)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added MPDBR on particles", len(plist1), ":", len(plist2), "at x0:kappa", str(x0), ":",
              str(kappa))

# wrapper for the SAMGR defined in custom c++ code
class SingleAxisMinGaussianRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, mean, sigma, axis, label, weight=1):
        particles = plist
        name = 'SingleAxisMinGaussianRestraint%1%'
        super(SingleAxisMinGaussianRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.desmosome.SingleAxisMinGaussianRestraint(particles, axis, mean, sigma)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added SAMGR on particles", len(plist), "at mean:sigma:axis", str(mean), ":", str(sigma), ":",
              str(axis))

# wrapper for the SPMGR defined in custom c++ code
class SinglePointMinGaussianRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, center, sigma, label, weight=1):
        particles = plist
        name = 'SinglePointMinGaussianRestraint%1%'
        super(SinglePointMinGaussianRestraintC, self).__init__(model, name=name, label=label, weight=weight)
        res_main = IMP.centrioles.SinglePointMinGaussianRestraint(particles, center, sigma)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added SPMGR on particles", len(plist), "at center:sigma", str(center), ":", str(sigma))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=data_directory+"PDB/",
                                  fasta_dir=data_directory+"FASTA/",
                                  gmm_dir=data_directory+"gmm/")

# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans = 0.55,
                                  max_rb_rot = 0.3,
                                  max_bead_trans = 4.25,
                                  max_srb_trans= 0.3,
                                  max_srb_rot=0.1)

# It's useful to have a list of the molecules.
molecules = t.get_components()

##### Uncomment the following lines to get test.rmf file to visualise the system representation
#
# # Uncomment this line for verbose output of the representation
IMP.atom.show_with_representations(root_hier)
# # output to RMF
fname = 'test_'+str(runID)+'.rmf'
rh = RMF.create_rmf_file(fname)
IMP.rmf.add_hierarchy(rh, root_hier)
IMP.rmf.save_frame(rh)


#####################################################
##################### RESTRAINTS ####################
#####################################################
# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object
# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of all molecule chains to included_objects to apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)
print("Excluded volume restraint applied")


# -----------------------------
# %%%%% MINIMUM PAIR DISTANCE BINDING RESTRAINT (MPDBR)
# MPDBR on Bld10 (1-652) - Sas4 (513-946)
for i in [0,1]:
    selection_tuple = (1, 652, 'Bld10', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    selection_tuple = (513, 946, 'Sas4', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    mpdbr1 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "BLD10-SAS4",
                                                              MPDBR_WEIGHT * Y2H_COIP_hetero_FACTOR)
    output_objects.append(mpdbr1)
print("MPDBR on Bld10 (1-652) to any Sas4 (513-946) applied")

for i in [0,1]:
    selection_tuple = (1, 652, 'Bld10', i, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    selection_tuple = (513, 946, 'Sas4', None, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    mpdbr2 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "BLD10-SAS4",
                                                              MPDBR_WEIGHT * Y2H_COIP_hetero_FACTOR)
    output_objects.append(mpdbr2)
print("MPDBR on Sas4 (513-946) to any Bld10 (1-652) applied")


# MPDBR on Bld10 (1-652) - Bld10 (1-652) applied
selection_tuple = (1, 652, 'Bld10', 0, None)
plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
selection_tuple = (1, 652, 'Bld10', 1, None)
plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
mpdbr3 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "BLD10-BLD10",
                                                      MPDBR_WEIGHT * Y2H_COIP_homo_FACTOR)
output_objects.append(mpdbr3)
print("MPDBR on Bld10 (1-652) to any Bld10 (1-652) applied")


# MPDBR on Bld10 (653-1298) - Bld10 (1-1298) applied
for i in [0,1]:
    selection_tuple = (653, 1298, 'Bld10', None, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 1298, 'Bld10', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    mpdbr4 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "BLD10-BLD10",
                                                      MPDBR_WEIGHT * Y2H_COIP_homo_FACTOR)
    output_objects.append(mpdbr4)
print("MPDBR on Bld10 (653-1298) to any Bld10 (1-1298) applied")

for i in [0,1]:
    selection_tuple = (653, 1298, 'Bld10', i, None)
    plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    selection_tuple = (1, 1298, 'Bld10', None, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
    mpdbr5 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "BLD10-BLD10",
                                                      MPDBR_WEIGHT * Y2H_COIP_homo_FACTOR)
    output_objects.append(mpdbr5)
print("MPDBR on Bld10 (1-1298) to any Bld10 (653-1298) applied")


# Sas4 (513-946) - Sas4 (513-946)
selection_tuple = (513, 946, 'Sas4', 0, None)
plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
selection_tuple = (513, 946, 'Sas4', 1, None)
plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
mpdbr6 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "SAS4-SAS4",
                                                      MPDBR_WEIGHT * Y2H_COIP_homo_FACTOR)
output_objects.append(mpdbr6)
print("MPDBR on Sas4 (513-946) to any Sas4 (513-946) applied")


# Sas4 (948-1342) - Sas4 (948-1342)
selection_tuple = (948, 1342, 'Sas4', 0, None)
plist1 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
selection_tuple = (948, 1342, 'Sas4', 1, None)
plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, MPDBR_RES)
mpdbr7 = MinimumPairDistanceBindingRestraint(mdl, plist1, plist2, 0, 1, "SAS4-SAS4",
                                                      MPDBR_WEIGHT * Y2H_COIP_homo_FACTOR)
output_objects.append(mpdbr7)
print("MPDBR on Sas4 (948-1342) to any Sas4 (948-1342) applied")


# # -------------------------
 # %%%%% SINGLE POINT MINIMUM GAUSSIAN RESTRAINT RESTRAINT (SPMGR)
 #
for i in [0,1]:
    selection_tuple = (1, 180, 'Bld10', i, None)
    plist = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple, SPMGR_RES)
    selection_tuple2 = (513, 664, 'Sas4', i, None)
    plist2 = IMP.pmi.tools.select_by_tuple_2(root_hier, selection_tuple2, SPMGR_RES)
    spmgr1 = SinglePointMinGaussianRestraintC(mdl, plist, [325, 250, 310], 10, 'SPMGRBld10-' + str(i), SPMGR_WEIGHT)
    output_objects.append(spmgr1)
    spmgr2 = SinglePointMinGaussianRestraintC(mdl, plist, [325, 250, 240], 10, 'SPMGRSas4-' + str(i), SPMGR_WEIGHT)
    output_objects.append(spmgr2)


print("Applied SPMGR on Bld10 and Sas4")


# # -------------------------
 # %%%%% EM RESTRAINT
 #
 # Scores a model based on its cross-correlation to an EM density.
 # Since cross-sorrelation is very expensive, we approximate both
 # the EM map and model as a set of 3D Gaussians (done in Representation).
 #
 # First, collect all density particles from the model.
 #Select residudes from Bld10.1 and Bld10.2; Sas4.1 and Sas4.2 to apply the restraint on
densities_A = IMP.atom.Selection(root_hier,molecule = "Bld10", residue_indexes = range(1,180), copy_indexes = [0,1], representation_type=IMP.atom.DENSITIES)
densities_B = IMP.atom.Selection(root_hier,molecule = "Sas4", residue_indexes = range(513,664), copy_indexes = [0,1], representation_type=IMP.atom.DENSITIES)
densities_AB = (densities_A | densities_B).get_selected_particles()

emr1 = IMP.pmi.restraints.em.GaussianEMRestraint(
             densities_AB,                 # Evaluate the restraint using these model densities
             target_fn=gmm_data_D,        # The EM map, approximated as a gaussian mixture model (GMM)
             slope=0.0000001,          # a small linear restraint to pull objects towards the EM map center
             weight=EM_WEIGHT)
output_objects.append(emr1)
print("EM Restraint Applied on Bld10 and Sas4")


# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are connected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()       #Add restraint to model
    output_objects.append(cr)

print("Connectivity restraint applied")


#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
#print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))
# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

#Shuffle MT binding regions at lower max_translation
MT_max_trans = 40
nonMT_max_trans = 60
#Select the MT binding residues
Sas4_MT_res = IMP.pmi.tools.select_by_tuple_2(root_hier, (513, 664, 'Sas4', None, None), 10)
Bld10_MT_res = IMP.pmi.tools.select_by_tuple_2(root_hier, (1, 180, 'Bld10', None, None), 10)

#Selecting nonMT binding residues
Sas4_nonMT_resN = IMP.pmi.tools.select_by_tuple_2(root_hier, (1, 512, 'Sas4', None, None), 30)
Sas4_nonMT_resC = IMP.pmi.tools.select_by_tuple_2(root_hier, (665, 1618, 'Sas4', None, None), 30)
Bld10_nonMT_res = IMP.pmi.tools.select_by_tuple_2(root_hier, (181, 1640, 'Bld10', None, None), 30)


MT_binding_res = [Bld10_MT_res, Sas4_MT_res]
print(MT_binding_res)
nonMT_binding_res = [Bld10_nonMT_res, Sas4_nonMT_resN, Sas4_nonMT_resC]

IMP.pmi.tools.shuffle_configuration(MT_binding_res,
                                    max_translation = MT_max_trans,
                                    excluded_rigid_bodies = nonMT_binding_res,             #Don't shuffle these rigid body objects
                                    hierarchies_included_in_collision = nonMT_binding_res) #Hierarchies that are not shuffled, but should be included in collision calculation (for fixed regions)

print("Shuffled beads in MT-binding region at max_trans = ",MT_max_trans)

#Shuffle non-MT binding regions at higher max_translation
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation = nonMT_max_trans,
                                    excluded_rigid_bodies = MT_binding_res,             #Don't shuffle these rigid body objects
                                    hierarchies_included_in_collision = MT_binding_res) #Hierarchies that are not shuffled, but should be included in collision calculation (for fixed regions)
print("Shuffled beads in non-MT-binding region at max_trans = ",nonMT_max_trans)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)
print(output_objects)
for i in output_objects:         #Add all restarints to the model
    i.add_to_model()

#print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = max_temp,
	    monte_carlo_sample_objects=dof.get_movers(),      #Pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory=run_output_dir,           #The output directory for this sampling run.
        output_objects=output_objects,                    #Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                             #Number of MC steps between writing frames
        number_of_best_scoring_models=0,                  #set >0 to store best PDB files (but this is slow)
        number_of_frames=num_frames)                      #Total number of frames to run / write to the RMF file.
        #test_mode=test_mode)                             #(Ignore this) Run in test mode (don't write anything)

# Ok, now we finally do the sampling!
rex.execute_macro()


# Outputs are then analyzed in a separate analysis script.
