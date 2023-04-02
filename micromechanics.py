#micromechanics module for myComposite Package create for ME227 class
#This module cover calculation of a lamina properties based on Rule of mixture
#Assume fiber and matrix are isotropic material
#Unidirectional Continuous fiber lamina only (Rule of Mixture)

import numpy as np

def calculate_E1(E_f:float, E_m:float ,V_f:float):
    #This function calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation
    #Inputs
    #E_f: Young's Modulus of the fiber
    #E_m: Young's Modulus of the matrix
    #V_f: Volume fraction of fiber
    return E_f*V_f + E_m(1-V_f)

def calculate_E2(E_f:float, E_m:float ,V_f:float):
    #This function calculate modulus in the traverse direction
    #Inputs
    #E_f: Young's Modulus of the fiber
    #E_m: Young's Modulus of the matrix
    #V_f: Volume fraction of fiber
    return (E_f*E_m)/(E_m*V_f + (1-V_f)*E_f)

def calculate_G(E:float, nu:float):
    #This function calculate shear modulus of an isotropic material
    #Inputs
    #E: Young's Modulus of the material
    #nu: Poisson's ratio of the material
    return E/(2*(1-nu))

def calculate_G(G_f:float, G_m:float, V_f:float):
    #This function calculate shear modulus of the lamaina
    #Inputs
    #G_f: Shear Modulus of the fiber
    #G_m: Shear Modulus of the matrix
    #V_f: Volume fraction of fiber
    return (G_f*G_m)/(G_m*V_f + (1-V_f)*G_f)

class Material:
    def __init__(self,name: str, modulus: float, poisson_ratio:float) -> None:
        self.name = name
        self.modulus = modulus
        self.poisson_ratio = poisson_ratio
        pass

class lamina:
    def __init__(self,fiberMaterial,) -> None:
        pass

