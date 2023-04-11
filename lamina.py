#Lamina module for myComposite Package create for ME227 class
#This module cover micromechanics and macromechanics calculation of a lamina
#Assume fiber and matrix are isotropic material
#Unidirectional Continuous fiber lamina only (Rule of Mixture)

import numpy as np
from dataclasses import dataclass, field

# Recommended SI units for in puts
# Modulus : Pa
# angle : deg
# Strength : pa
# thickness : m
# N(force per unit width) = N/m
# M(moment per init width) = N m / m

def calculate_E1(E_f:float, E_m:float ,V_f:float) -> float: 
    """Calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation based on Rule of mixture
    Inputs
    E_f: Young's Modulus of the fiber
    E_m: Young's Modulus of the matrix
    V_f: Volume fraction of fiber"""
    return E_f*V_f + E_m*(1-V_f)

def calculate_E2(E_f:float, E_m:float ,V_f:float) -> float:
    """Calculate modulus in the traverse direction
    Inputs
    E_f: Young's Modulus of the fiber
    E_m: Young's Modulus of the matrix
    V_f: Volume fraction of fiber"""
    return (E_f * E_m)/(E_m * V_f + (1-V_f) * E_f)

def calculate_nu12 (nu_f:float, nu_m:float ,V_f:float) -> float:
    """Calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation based on Rule of mixture
    Inputs
    nu_f: Poisson's ratio of the fiber
    nu_m: Poisson's ratio of the matrix
    V_f: Volume fraction of fiber"""
    return nu_f*V_f + nu_m*(1-V_f)

def calculate_material_G(E:float, nu:float) -> float:
    """Calculate shear modulus of an isotropic material 
    Inputs
    E: Young's Modulus of the material
    nu: Poisson's ratio of the material"""
    return E/(2*(1-nu))

def calculate_material_nu21(nu_12:float,E_1:float,E_2:float)-> float:
    """Calculate shear modulus of an isotropic material 
    Inputs
    nu_12: Poisson's ratio of the material
    E_1: Young's Modulus of the material in  dir 1
    E_2: Young's Modulus of the material in  dir 2"""
    return nu_12*E_2/E_1

def calculate_G12(G_f:float, G_m:float, V_f:float) -> float:
    """Calculate shear modulus of the lamaina based on Rule of mixture
    Inputs
    G_f: Shear Modulus of the fiber
    G_m: Shear Modulus of the matrix
    V_f: Volume fraction of fiber"""
    return (G_f*G_m)/(G_m*V_f + (1-V_f)*G_f)

#TODO: Halpin-Tsai Model for UD and random and Randomly Orient
#Halpin-Tsai : Unidirectional Continuous
def calculate_halpinTsai(E_m,E_f,V_f,psi):
  E_ratio = E_f/E_m
  eta = (E_ratio-1)/(E_ratio+psi)
  return (E_m*(1+psi*eta*V_f))/(1-eta*V_f)

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

def assemble_transformation_matrix(angle):
    """ Assembles T matrix (angles transformation matrix LCS -> MCS).
    Inputs
    angle: rotation angle (deg)"""
    if not isinstance(angle, (int, float)) or not (-360 <= angle <= 360):
        raise Exception("lamina angle is not between +- 360 degrees"+
                                " or it is not an int/float")

    #Transforms angle (degrees) to angle (radians)
    angle = np.pi*angle/180

    cos = np.cos(angle)
    sin = np.sin(angle)
    cs = cos*sin
    cc = cos**2
    ss = sin**2

    T = np.zeros((3, 3))
    T = np.array([  [cc,    ss,   2*cs   ],
                    [ss,    cc,  -2*cs   ],
                    [-cs, cs, cc-ss]])
    return T

def assemble_reduced_stiffness_matrix (E_1,E_2, G_12, nu_12, nu_21: float = None):
    """ Assembles Q matrix (reduced stiffness matrix) for a lamina. 
    Inputs
    E_1: Young's Modulus of the lamina
    E_2: Young's Modulus of the lamina
    nu_12: Poisson's ratio of the lamina
    G_12: Shear Modulus of the lamina"""
    if nu_21 == None:
        nu_21 = calculate_material_nu21(nu_12,E_1,E_2)

    Q11 = E_1/(1-nu_12*nu_21)        
    Q12 = nu_12*E_1*E_2/(E_1-(nu_12**2)*E_2)
    Q22 = E_1*E_2/(E_1-(nu_12**2)*E_2)
    Q66 = G_12

    Q = np.zeros((3, 3))
    Q = np.array([[Q11, Q12, 0],
                     [Q12, Q22, 0],
                     [0,   0, Q66]])
    return Q

def calculate_matrixQbar (matrixQ, angle, matrixT = None):
    """ Calculate Qbar matrix (transformed reduced stiffness matrix) for a lamina. 
    Inputs
    matrixQ: reduced stiffness matrix in principal fiber direction
    angle: angle: rotation angle (deg)
    matrixT: angles transformation matrix. If defined, the T matrix calculation will be skipped.  """

    if matrixT is not None:
        if not isinstance(matrixT, (np.ndarray)):
            raise Exception("angles transformation matrix need to be a 3x3 array")
    else:
        T = assemble_transformation_matrix(angle)
    T_inv = np.linalg.inv(T)
    Qbar = T_inv@ matrixQ @ T_inv.T
    return Qbar

@dataclass
class Strength:
    """Class for containing description of a laminate strength."""
    S_Lt:float
    S_Lc:float
    S_Tt:float
    S_Tc:float
    S_shear: float

class Material:
    def __init__(self, E_1: float, E_2: float, nu_12, G_12: float, strength:Strength = None , name:str = "undefined", is_isotropic: bool = False) -> None:
        self.name = name
        self.E_1 = float(E_1)
        self.E_2 = float(E_2)
        self.nu_12 = float(nu_12)
        self.nu_21 = calculate_material_nu21(self.nu_12,self.E_1,self.E_2)
        self.G_12 = G_12
        self.strength  = strength
        self.matrixQ = assemble_reduced_stiffness_matrix(self.E_1,self.E_2, self.G_12, self.nu_12,self.nu_21)
        self.is_isotropic = is_isotropic
        
    @classmethod
    def isotropic(cls,E: float,nu: float,G: float = None, strength:Strength = None ,name:str = "undefined") -> None:
        if G is None:
            G = calculate_material_G(E,nu)
        return cls(E_1= E, E_2= E, nu_12= nu, G_12 =G, strength = strength, name= name, is_isotropic = True)
    
    def __str__(self) -> str:
        return "material name: %s, E_1: %s, E_2: %s, nu_12: %s, nu_21: %s, G_12: %s, is_isotropic: %s" % (self.name, self.E_1, self.E_2, self.nu_12, self.nu_21, self.G_12, self.is_isotropic)
    
class Lamina(Material):    
    def __init__(self,fiber: Material, matrix: Material, V_f: float, thickness:float, strength:Strength = None , density:float = None, name:str = "undefined") -> None:
        if not (fiber.is_isotropic and matrix.is_isotropic):
            raise Exception("Fiber and Matrix must be isotropic material")
        #remember the fiber and matrix
        self.name = name
        self.fiber = fiber
        self.matrix = matrix
        
        self.E_1 = calculate_E1(fiber.E_1,matrix.E_1,V_f)
        self.E_2 = calculate_halpinTsai(fiber.E_1,matrix.E_1,V_f,2)
        self.nu_12 = calculate_nu12(fiber.nu_12,matrix.nu_12,V_f)
        self.nu_21 = calculate_material_nu21(self.nu_12,self.E_1,self.E_2)
        self.G_12 = calculate_G12(fiber.G_12,matrix.G_12,V_f)
        self.matrixQ = assemble_reduced_stiffness_matrix(self.E_1,self.E_2, self.G_12,self.nu_12)
        self.thickness = thickness
        self.strength  = strength
        self.density  = density
        self.is_isotropic = False

    @classmethod
    def from_Mat_props(cls, E_1: float, E_2: float, nu_12, G_12: float, thickness:float, strength:Strength = None, density:float = None,name:str = "undefined"):
        obj = cls.__new__(cls)
        super(Lamina, obj).__init__(E_1,E_2,nu_12,G_12,name)
        
        #remember the fiber and matrix
        obj.name = name
        obj.fiber = None
        obj.matrix = None

        obj.E_1 = float(E_1)
        obj.E_2 = float(E_2)
        obj.nu_12 = float(nu_12)
        obj.nu_21 = calculate_material_nu21(obj.nu_12,obj.E_1,obj.E_2)
        obj.G_12 = G_12
        obj.matrixQ = assemble_reduced_stiffness_matrix(obj.E_1,obj.E_2, obj.G_12, obj.nu_12,obj.nu_21)
        obj.thickness = thickness
        obj.strength  = strength
        obj.density  = density
        obj.is_isotropic = False
        return obj

    

