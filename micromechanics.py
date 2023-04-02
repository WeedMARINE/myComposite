#micromechanics module for myComposite Package create for ME227 class
#This module cover calculation of a lamina properties based on Rule of mixture
#Assume fiber and matrix are isotropic material
#Unidirectional Continuous fiber lamina only (Rule of Mixture)

import numpy as np

def calculate_E1(E_f:float, E_m:float ,V_f:float) -> float: 
    """Calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation
    Inputs
    E_f: Young's Modulus of the fiber
    E_m: Young's Modulus of the matrix
    V_f: Volume fraction of fiber"""
    return E_f*V_f + E_m(1-V_f)

def calculate_E2(E_f:float, E_m:float ,V_f:float) -> float:
    """Calculate modulus in the traverse direction
    Inputs
    E_f: Young's Modulus of the fiber
    E_m: Young's Modulus of the matrix
    V_f: Volume fraction of fiber"""
    return (E_f*E_m)/(E_m*V_f + (1-V_f)*E_f)

def calculate_nu12 (nu_f:float, nu_m:float ,V_f:float) -> float:
    """Calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation
    Inputs
    nu_f: Poisson's ratio of the fiber
    nu_m: Poisson's ratio of the matrix
    V_f: Volume fraction of fiber"""
    return nu_f*V_f + nu_m(1-V_f)

def calculate_G(E:float, nu:float) -> float:
    """Calculate shear modulus of an isotropic material
    Inputs
    E: Young's Modulus of the material
    nu: Poisson's ratio of the material"""
    return E/(2*(1-nu))

def calculate_G(G_f:float, G_m:float, V_f:float) -> float:
    """Calculate shear modulus of the lamaina
    Inputs
    G_f: Shear Modulus of the fiber
    G_m: Shear Modulus of the matrix
    V_f: Volume fraction of fiber"""
    return (G_f*G_m)/(G_m*V_f + (1-V_f)*G_f)



#TODO: Halpin-Tsai Model for UD and random and Randomly Orient
# #Halpin-Tsai : Unidirectional Continuous
# def cal_halpinTsai(Em,Ef,Vf,psi):
#   E_ratio = Ef/Em
#   eta = (E_ratio-1)/(E_ratio+psi)
#   return (Em*(1+psi*eta*Vf))/(1-eta*Vf)

# #Halpin-Tsai Unidirectional Calculation
# HT_E1 = cal_halpinTsai(3.5,72,0.4,100)
# HT_E2 = cal_halpinTsai(3.5,72,0.4,2)
# print(HT_E1)
# print(HT_E2)

# #Unidirectional Discontinuous
# def cal_MatProp_UD(Em,Ef,Vf,AR):
#   #This function calculate E1,E2 of a Unidirectional Discontinuous fiber lamina. It assume isotropy in Matrix and Fiber.
#   #TODO: add G12,v12 calculation
#   E_ratio = Ef/Em
#   etaL = (E_ratio-1)/(E_ratio+2*AR)
#   etaT = (E_ratio-1)/(E_ratio+2)
#   E1 = (Em*(1+2*AR*etaL*Vf))/(1-etaL*Vf)
#   E2 = (Em*(1+2*etaT*Vf))/(1-etaT*Vf)
#   return E1,E2

# #Randomly Oriented Discontinuous Fiber Lamina
# def cal_MatProp_RD(Em,Ef,Vf,AR):
#     #This function calculate E1,E2 of a Randomly Oriented Discontinuous Fiber Lamina. It assume isotropy in Matrix and Fiber.
#     #A thin lamina containing randomly oriented discontinuous fibers exhibits planar isotropic behavior
#     E1,E2 = cal_MatProp_UD(Em,Ef,Vf,AR)
#     E_rand = (3/8)*E1 + (5/8)*E2
#     G_rand = (1/8)*E1 + (1/4)*E2
#     v_rand = (E_rand/(2*G_rand))-1
#     return E_rand,G_rand,v_rand