from laminate import *

def maximum_stress_2D(lamina:Lamina, angle, stress_xy):
    """check if a lamina fail according to maximum_stress criteria"""
    is_failed = True
    stress_12 = assemble_transformation_matrix(angle) @ stress_xy
    
    working_stress = np.array(stress_12)
    working_stress[2] = np.abs(working_stress[2]) #avoid negative value

    strength = np.array([lamina.strength.S_Lt, lamina.strength.S_Tt, lamina.strength.S_shear])
    if working_stress[0] < 0:
        strength[0] = lamina.strength.S_Lc
    if working_stress[1] < 0:
        strength[0] = lamina.strength.S_Tc

    factor_of_safety = strength/np.abs(working_stress)

    if np.all(factor_of_safety >= 1):
        is_failed = False

    return is_failed, np.min(factor_of_safety)

def tsaihill_2D(lamina:Lamina, angle, stress_xy):
    """check if a lamina fail according to Tsai-Hill criteria"""
    is_failed = True
    stress_12 = assemble_transformation_matrix(angle) @ stress_xy
    
    working_stress = np.array(stress_12)
    working_stress[2] = np.abs(working_stress[2]) #avoid negative value

    strength = np.array([lamina.strength.S_Lt, lamina.strength.S_Tt, lamina.strength.S_shear])
    if working_stress[0] < 0:
        strength[0] = lamina.strength.S_Lc
    if working_stress[1] < 0:
        strength[0] = lamina.strength.S_Tc
    
    failure_index = (working_stress[0]/strength[0])**2 - (working_stress[0]*working_stress[1]/strength[0]**2) + (working_stress[1]/strength[1])**2 + (working_stress[2]/strength[2])**2

    factor_of_safety = 1/np.sqrt(failure_index)
    if factor_of_safety >= 1:
        is_failed = False

    return is_failed, factor_of_safety

def laminate_failure_check(laminate:Laminate, failure_criteria: callable, load):
    """check if a laminate fail according to maximum_stress criteria"""
    matrix_abd = np.linalg.inv(laminate.matrixABD)
    vector_strain_curve = matrix_abd @ load
    midplane_strain_xy = vector_strain_curve[:3]
    midplane_curvature = vector_strain_curve[3:]

    layer_count = len(laminate.vectorZ) - 1
    plies_is_failed = []
    plies_factor_of_safety = []

    for k in range(layer_count):
        ply_strain_xy = midplane_strain_xy + laminate.vectorZ[k+1]*midplane_curvature
        ply_stress_xy = laminate.description.plies[k].matrixQ @ ply_strain_xy
        ply_stress_12 = assemble_transformation_matrix(-laminate.description.angles[k]) @ ply_stress_xy
        is_failed, factor_of_safety = failure_criteria(laminate.description.plies[k], laminate.description.angles[k],ply_stress_12)
        plies_is_failed.append(is_failed)
        plies_factor_of_safety.append(factor_of_safety)
    
    return plies_is_failed, plies_factor_of_safety
