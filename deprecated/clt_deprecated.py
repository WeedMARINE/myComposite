"""
lamipy project - laminated composites calculations in Python.

clt.py - Module containing functions for stress&strain calculations 
according to the CLASSICAL LAMINATE THEORY. 

INPUTS:
lam -- laminate layup characteristics
mat_list -- list with dictionaries of material properties
F -- 'Force' vector (generalized stress vector) contains applied loads on 
the composite
fail_list -- list with failed layers and their failure modes
dT -- temperature variation
dM -- moisture variation

MAIN FUNCTION:
calc_stressCLT -- Calculates stress and strain vectors according to the
Classical Laminate Theory, at the top and bottom of each layer.

OUTPUTS:
Dictionary with:
Material Coord. System stresses and strains, top and bottom of layer;
Laminate Coord. System strains, top and bottom of layer.

Joao Paulo Bernhardt - October 2017
"""

import numpy
import macromechanics

# class LaminateLayupError(TypeError): pass

# def calc_stressCLT(mat_list, lam, F, fail_list = None, dT = 0, dM = 0):
#     """ Calculates stress and strain vectors according to CLT. """

#     # Get number of layers
#     num = len(lam["ang"])
#     total_thk = numpy.sum(lam["thk"])

#     # Z vector contains z coordinates.
#     Z = assemble_Z(lam)

#     # Calculates stiffness matrix based on laminate properties.
#     ABD = assemble_ABD(mat_list, lam, Z, fail_list)

#     # Calculates thermal, moisture resultants and adds to Forces vector
#     Nt = calc_thermal_forces(mat_list, lam, Z, fail_list, dT)
#     Nm = calc_moisture_forces(mat_list, lam, Z, fail_list, dM)
#     F = F + Nt + Nm

#     # Calculates strain vector by solving eq. form AX = B
#     strain_vector = numpy.linalg.solve(ABD, F)

#     # Initializes strain and curvature vectors & defines values
#     strains = curvatures = numpy.zeros((1, 3))
#     strains = strain_vector[:3]
#     curvatures = strain_vector [3:6]

#     # Initializes Laminate System (LS) strain vectors (inferior and superior)
#     LS_strain_inf = numpy.zeros((3, num))
#     LS_strain_sup = numpy.zeros((3, num))

#     # Assign Laminate System strain (Epsilon = Epsilon0 + kappa*z)
#     for i in range(num):
#         LS_strain_inf[:3,i] = strains + curvatures*Z[i]
#         LS_strain_sup[:3,i] = strains + curvatures*Z[i+1]

#     # Initialize Material System stress & strain vectors (inferior and superior)
#     MS_strain_inf = numpy.zeros((3, num))
#     MS_strain_sup = numpy.zeros((3, num))
#     MS_stress_inf = numpy.zeros((3, num))
#     MS_stress_sup = numpy.zeros((3, num))

#     # Calculates Material System stresses and strains
#     for i in range(num):
#         mat_id = lam["mat_id"][i]
#         mat_prop = mat_list[mat_id]
        
#         # Checks if there's failure list
#         if isinstance(fail_list, list):
#             Q = assemble_matrixQ(mat_prop, fail_list[i])
#         else:
#             Q = assemble_matrixQ(mat_prop)
        
#         T =  macromechanics.assemble_matrixT(lam["ang"][i])
#         MS_strain_sup[:,i] = numpy.dot(T, LS_strain_sup[:,i])
#         MS_strain_inf[:,i] = numpy.dot(T, LS_strain_inf[:,i])
#         MS_stress_sup[:,i] = numpy.dot(Q, MS_strain_sup[:,i])
#         MS_stress_inf[:,i] = numpy.dot(Q, MS_strain_inf[:,i])

#     # Outputs the stresses and strains vectors.
#     # LCS = Laminate System; MCS = Material Coordinate System; 
#     # inf = inferior; sup = superior.
#     return {"LCS" : {"strain" : {"inf" : LS_strain_inf,
#                                  "sup" : LS_strain_sup}}, 
#             "MCS" : {"stress" : {"inf" : MS_stress_inf,
#                                  "sup" : MS_stress_sup}, 
#                      "strain" : {"inf" : MS_strain_inf,
#                                  "sup" : MS_strain_sup}}}


# ##### BELOW: Auxiliary functions for the calculation of stresses and strains

# def assemble_Z(lam):
#     """ Assembles Z vector which contains z coordinates of laminate. """

#     if not isinstance(lam, dict):
#         raise LaminateLayupError('Laminate has to be a dictionary.')

#     # Get number of layers
#     num = len(lam["ang"])
#     total_thk = numpy.sum(lam["thk"])
#     Z = numpy.zeros(num + 1)
#     Z[0] = -total_thk/2

#     for i in range(num):
#         Z[i + 1] = Z[i] + lam["thk"][i]

#     return Z

# def calc_thermal_forces(mat_list, lam, Z, fail_list = None, dT = 0):
#     """ Calculates force resultants due to temperature variations. """
    
#     if (not (isinstance(mat_list, (tuple, list))) or
#     not (isinstance(lam, dict)) or
#     not (isinstance(Z, numpy.ndarray))):
#         raise LaminateLayupError('invalid input(s) for this function')

#     num = len(lam["ang"])

#     Nt = numpy.zeros(6)
#     a = numpy.zeros(3)
#     a_LCS = numpy.zeros(3)
    

#     for i in range(num):
#         mat_id = lam["mat_id"][i]
#         mat_prop = mat_list[mat_id]
#         ang = lam["ang"][i]        

#         #alpha (material coord system) vector
#         a = [mat_prop["a1"], mat_prop["a2"], 0]

#         T =  macromechanics.assemble_matrixT(ang)
#         Ti = numpy.transpose(T)

#         a_LCS = numpy.matmul(Ti, a)

#         # Checks if there's failure list
#         if isinstance(fail_list, list):
#             Q = assemble_matrixQ(mat_prop, fail_list[i])
#         else:
#             Q = assemble_matrixQ(mat_prop)

#         Qbar = numpy.dot(numpy.dot(Ti, Q), T)

#         Nt[:3] = Nt[:3] + numpy.dot(Qbar, a_LCS) * dT * (Z[i+1] - Z[i])
#         Nt[3:6] = Nt[3:6] + numpy.dot(Qbar, a_LCS) * dT * (1/2) * \
#                                                      (Z[i+1]**2 - Z[i]**2)

#     return Nt 

# def calc_moisture_forces(mat_list, lam, Z, fail_list = None, dM = 0):
#     """ Calculates force resultants due to moisture variations. """

#     if (not (isinstance(mat_list, (tuple, list))) or
#     not (isinstance(lam, dict)) or
#     not (isinstance(Z, numpy.ndarray))):
#         raise LaminateLayupError('invalid input(s) for this function')


#     num = len(lam["ang"])

#     Nm = numpy.zeros(6)
#     b = numpy.zeros(3)
#     b_LCS = numpy.zeros(3)
    

#     for i in range(num):
#         mat_id = lam["mat_id"][i]
#         mat_prop = mat_list[mat_id]
#         ang = lam["ang"][i]        

#         #alpha (material coord system) vector
#         b = [mat_prop["b1"], mat_prop["b2"], 0]

#         T = macromechanics.assemble_matrixT(ang)
#         Ti = numpy.transpose(T)

#         b_LCS = numpy.matmul(Ti, b)

#         # Checks if there's failure list
#         if isinstance(fail_list, list):
#             Q = assemble_matrixQ(mat_prop, fail_list[i])
#         else:
#             Q = assemble_matrixQ(mat_prop)

#         Qbar = numpy.dot(numpy.dot(Ti, Q), T)

#         Nm[:3] = Nm[:3] + numpy.dot(Qbar, b_LCS) * dM * (Z[i+1] - Z[i])
#         Nm[3:6] = Nm[3:6] + numpy.dot(Qbar, b_LCS) * dM * (1/2) * \
#                                                 (Z[i+1]**2 - Z[i]**2)

#     return Nm 

# def assemble_ABD(mat_list, lam, Z, fail_list = None):
#     """ Assembles ABD matrix (laminate stiffness matrix). """

#     if (not (isinstance(mat_list, (tuple, list))) or
#     not (isinstance(lam, dict)) or
#     not (isinstance(Z, numpy.ndarray))):
#         raise LaminateLayupError('invalid input(s) for this function')

#     num = len(lam["ang"])
#     A = B = D = numpy.zeros((3,3))

#     for i in range(num):
#         mat_id = lam["mat_id"][i]
#         mat_prop = mat_list[mat_id]
#         ang = lam["ang"][i]        

#         T = macromechanics.assemble_matrixT(ang)
#         Ti = numpy.transpose(T)

#         # Checks if there's failure list
#         if isinstance(fail_list, list):
#             Q = assemble_matrixQ(mat_prop, fail_list[i])
#         else:
#             Q = assemble_matrixQ(mat_prop)

#         Qi = numpy.matmul(numpy.matmul(Ti, Q), T)

#         A = A + Qi*(Z[i+1] - Z[i])
#         B = B + (1/2)*Qi*(Z[i+1]**2 - Z[i]**2)
#         D = D + (1/3)*Qi*(Z[i+1]**3 - Z[i]**3)

#     ABD = numpy.zeros((6,6))
#     ABD[:3, :3] = A
#     ABD[:3, 3:6] = ABD[3:6, :3] = B
#     ABD[3:6, 3:6] = D

#     return ABD