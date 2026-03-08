from openawsem.functionTerms import *
from openawsem.helperFunctions.myFunctions import *

try:
    from openmm.unit import angstrom
    from openmm.unit import kilocalorie_per_mole
except ModuleNotFoundError:
    from simtk.unit import angstrom
    from simtk.unit import kilocalorie_per_mole

def set_up_forces(oa, computeQ=False, submode=-1, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = [
        basicTerms.con_term(oa),
        basicTerms.chain_term(oa),
        basicTerms.chi_term(oa),
        basicTerms.excl_term(oa),
        basicTerms.rama_term(oa),
        basicTerms.rama_proline_term(oa),
        basicTerms.rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        contactTerms.contact_term(oa, k_contact=1*4.184, k_burial=1*4.184),
        hydrogenBondTerms.beta_term_1(oa, version='efficiency_optimized', k=0.5*4.184),
        hydrogenBondTerms.beta_term_2(oa, version='efficiency_optimized', k=0.5*4.184),
        hydrogenBondTerms.beta_term_3(oa, version='efficiency_optimized', k=0.5*4.184),
        hydrogenBondTerms.pap_term_1(oa, version='efficiency_optimized', k=0.5*4.184),
        hydrogenBondTerms.pap_term_2(oa, version='efficiency_optimized', k=0.5*4.184),
        # templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=False),
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False,
                                               k_fm=0.01*4.184, fm_well_width=0.1, frag_table_dr=0.01, min_seq_sep=3, max_seq_sep=9),
        debyeHuckelTerms.debye_huckel_term(oa, chargeFile="charge.txt"),
    ]
    if computeQ:
        forces.append(biasTerms.rg_term(oa))
        forces.append(biasTerms.q_value(oa, "crystal_structure-cleaned.pdb", forceGroup=1))
        forces.append(biasTerms.qc_value(oa, "crystal_structure-cleaned.pdb"))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", residueIndexGroup=list(range(0, 15)), forceGroup=1))
    if submode == 0:
        additional_forces = [
            # contact_term(oa),
        ]
        forces += additional_forces
    return forces
