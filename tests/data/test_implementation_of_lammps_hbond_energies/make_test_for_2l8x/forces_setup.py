from openawsem.functionTerms import *
from openawsem.helperFunctions.myFunctions import *

try:
    from openmm.unit import angstrom
    from openmm.unit import kilocalorie_per_mole
except ModuleNotFoundError:
    from simtk.unit import angstrom
    from simtk.unit import kilocalorie_per_mole

def set_up_forces(oa, computeQ=False, submode=-1, contactParameterLocation=".", membrane_center=-0*angstrom, extension="pdb"):
    # apply forces
    forces = [
        basicTerms.con_term(oa),
        basicTerms.chain_term(oa),
        basicTerms.chi_term(oa),
        basicTerms.excl_term(oa, periodic=False),
        basicTerms.rama_term(oa),
        basicTerms.rama_proline_term(oa,k_rama_proline=3*4.184),
        basicTerms.rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        contactTerms.contact_term(oa,k_contact=0.75*4.184,k_burial=4.184),
        #hydrogenBondTerms.beta_term_1(oa,forceGroup=23,k=0.5*4.184),
        #hydrogenBondTerms.beta_term_2(oa,forceGroup=24,k=0.5*4.184),
        #hydrogenBondTerms.beta_term_3(oa,forceGroup=25,k=0.5*4.184),
        hydrogenBondTerms.beta_term_1_old(oa,forceGroup=23,k_beta=0.5*4.184),
        hydrogenBondTerms.beta_term_2_old(oa,forceGroup=24,k_beta=0.5*4.184),
        hydrogenBondTerms.beta_term_3_old(oa,forceGroup=25,k_beta=0.5*4.184),
        #hydrogenBondTerms.pap_term_1(oa),
        #hydrogenBondTerms.pap_term_2(oa),
        hydrogenBondTerms.pap_term_old(oa),
        # membraneTerms.membrane_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center),
        # membraneTerms.membrane_preassigned_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center, zimFile="PredictedZim"),
        # templateTerms.er_term(oa),
        # templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        #templateTerms.fragment_memory_term(oa, frag_file_list_file="./fragsLAMW.mem.single", npy_frag_table="./single_frags.npy", UseSavedFragTable=False,k_fm=20*0.04184,fm_well_width=0.01,forceGroup=29),
        debyeHuckelTerms.debye_huckel_term(oa, chargeFile="charge.txt"),
    ]
    if computeQ:
        forces.append(biasTerms.rg_term(oa))
        forces.append(biasTerms.q_value(oa, f"crystal_structure-cleaned.{extension}", forceGroup=1))
        forces.append(biasTerms.qc_value(oa, f"crystal_structure-cleaned.{extension}"))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", residueIndexGroup=list(range(0, 15)), forceGroup=1))
    if submode == 0:
        additional_forces = [
            # contact_term(oa),
        ]
        forces += additional_forces
    return forces
