try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
import numpy as np

def membrane_term(oa, k=1*kilocalorie_per_mole, k_m=20, z_m=1.5, membrane_center=0*angstrom, forceGroup=24):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_membrane = k * oa.k_awsem

    # 02082024 Rebekah Added --- End
    if oa.periodic_box:
        membrane = CustomExternalForce(f"k_membrane * hydrophobicityScale *"
                                       f"( 0.5*tanh({k_m}*(z_periodic + {z_m})) + 0.5*tanh({k_m}*({z_m} - z_periodic)) )"
                                       f";z_periodic=periodicdistance(0,0,z,0,0,{membrane_center})")
        print("Membrane_term is Periodic")
    else:
        membrane = CustomExternalForce (f"k_membrane*\
                (0.5*tanh({k_m}*((z-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z-{membrane_center}))))*hydrophobicityScale")                          
    # 02082024 Rebekah Added --- End

    membrane.addPerParticleParameter("hydrophobicityScale")
    membrane.addGlobalParameter("k_membrane", k_membrane)
    zim = np.loadtxt("zim")
    
    # cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    ca = oa.ca
    for i in ca:
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
        membrane.addParticle(i, [zim[oa.resi[i]]])
    membrane.setForceGroup(forceGroup)
    return membrane



def membrane_with_pore_term(oa, k=1*kilocalorie_per_mole, pore_center_x=0, pore_center_y=0, pore_radius=10, k_pore=0.1, k_m=20, z_m=1.5, membrane_center=0*angstrom, forceGroup=24):
    # inside the pore, the energy is zero, as if the residue is in the water.
    # pore_center_x, pore_center_y, pore_radius in unit of nanometer.

    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_membrane = k * oa.k_awsem

    membrane = CustomExternalForce(f"{k_membrane}*\
            (0.5*tanh({k_m}*((z-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z-{membrane_center}))))*(1-alpha)*hydrophobicityScale;\
            alpha=0.5*(1+tanh({k_pore}*({pore_radius}-rho)));\
            rho=((x-{pore_center_x})^2+(y-{pore_center_y})^2)^0.5")

    membrane.addPerParticleParameter("hydrophobicityScale")
    zim = np.loadtxt("zim")
    # cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    ca = oa.ca
    for i in ca:
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
        membrane.addParticle(i, [zim[oa.resi[i]]])
    membrane.setForceGroup(forceGroup)
    return membrane


def membrane_preassigned_term(oa, k=1*kilocalorie_per_mole, k_m=20, z_m=1.5, membrane_center=0*angstrom, zimFile="zimPosition", forceGroup=24):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_membrane = k * oa.k_awsem

    membrane = CustomExternalForce(f"{k_membrane}*\
            (0.5*tanh({k_m}*((z-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z-{membrane_center}))))*zim")
    membrane.addPerParticleParameter("zim")
    # zim = np.loadtxt("zim")
    zimPosition = np.loadtxt(zimFile)
    zim = [-1 if z == 2 else 1 for z in zimPosition]
    # print(zim)
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    for i in cb_fixed:
        membrane.addParticle(i, [zim[oa.resi[i]]])
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    membrane.setForceGroup(forceGroup)
    return membrane


def SideToZ_m(side):
    side = side.strip()
    if side == "down":
        return -1.5
    if side == "up":
        return 1.5
    if side == "middle":
        return 0

def membrane_preassigned_side_term(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, zimFile="PredictedZimSide", forceGroup=24):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_membrane = k * oa.k_awsem
    membrane = CustomExternalForce(f"{k_membrane}*(abs(z-{membrane_center}-z_m))")
    membrane.addPerParticleParameter("z_m")

    with open(zimFile) as f:
        a = f.readlines()

    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    # print(cb_fixed)
    for i in cb_fixed:
        z_m = SideToZ_m(a[oa.resi[i]])
        # print(oa.resi[i])
        membrane.addParticle(i, [z_m])
    membrane.setForceGroup(forceGroup)
    return membrane


def single_helix_orientation_bias_term(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, z_m=1.5, k_m=20, atomGroup=-1, forceGroup=18):
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_single_helix_orientation_bias = oa.k_awsem * k
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup
    n = len(group)
    theta_z1 = f"(0.5*tanh({k_m}*((z1-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z1-{membrane_center}))))"
    theta_z2 = f"(0.5*tanh({k_m}*((z2-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z2-{membrane_center}))))"
    normalization = n * (n - 1) / 2
    v_orientation = CustomCompoundBondForce(2, f"helix_orientation*{k_single_helix_orientation_bias}/{normalization}*((x1-x2)^2+(y1-y2)^2)*{theta_z1}*{theta_z2}")
    v_orientation.addGlobalParameter("helix_orientation", 1)
    # rcm_square = CustomCompoundBondForce(2, "1/normalization*(x1*x2)")
    # v_orientation.addGlobalParameter("k_single_helix_orientation_bias", k_single_helix_orientation_bias)
    # rg_square = CustomBondForce("1/normalization*(sqrt(x^2+y^2)-rcm))^2")
    # rg = CustomBondForce("1")
    # v_orientation.addGlobalParameter("normalization", n*n)
    for i in group:
        for j in group:
            if j <= i:
                continue
            v_orientation.addBond([ca[i], ca[j]], [])

    v_orientation.setForceGroup(forceGroup)
    return v_orientation






