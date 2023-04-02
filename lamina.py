#Lamina module for myComposite Package create for ME227 class
#This module cover micromechanics and macromechanics calculation of a lamina
#Assume fiber and matrix are isotropic material
#Unidirectional Continuous fiber lamina only (Rule of Mixture)

import numpy as np

class LaminaCalculationError(TypeError): pass

def calculate_E1(E_f:float, E_m:float ,V_f:float) -> float: 
    """Calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation based on Rule of mixture
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
    """Calculate modulus of a unidirectional continuous lamina in the direction parallel to the fiber orientation based on Rule of mixture
    Inputs
    nu_f: Poisson's ratio of the fiber
    nu_m: Poisson's ratio of the matrix
    V_f: Volume fraction of fiber"""
    return nu_f*V_f + nu_m(1-V_f)

def calculate_material_G(E:float, nu:float) -> float:
    """Calculate shear modulus of an isotropic material 
    Inputs
    E: Young's Modulus of the material
    nu: Poisson's ratio of the material"""
    return E/(2*(1-nu))

def calculate_G_12(G_f:float, G_m:float, V_f:float) -> float:
    """Calculate shear modulus of the lamaina based on Rule of mixture
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

def assemble_matrixT(angle):
    """ Assembles T matrix (angles transformation matrix LCS -> MCS).
    Inputs
    angle: rotation angle (deg)"""
    if not isinstance(angle, (int, float)) or not (-360 <= angle <= 360):
        raise LaminaCalculationError("lamina angle is not between +- 360 degrees"+
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

def assemble_matrixQ (E_1,E_2, nu_12, G_12):
    """ Assembles Q matrix (reduced stiffness matrix) for a lamina. 
    Inputs
    E_1: Young's Modulus of the lamina
    E_2: Young's Modulus of the lamina
    nu_12: Poisson's ratio of the lamina
    G_12: Shear Modulus of the lamina"""

    Q11 = 1/E_1        
    Q12 = -nu_12/E_1
    Q22 = 1/E_2
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
            raise LaminaCalculationError("angles transformation matrix need to be a 3x3 array")
        T = matrixT
    else:
        T = assemble_matrixT(angle)
    T_inv = np.linalg.inv(T)
    Qbar = T_inv@ matrixQ @ T_inv.T
    return Qbar

class Material:
    def __init__(self,E,nu,G: float = None) -> None:
        self.E = E
        self.nu = nu
        if G is None:
            self.G = calculate_material_G(E,nu)
        else:
            self.G = G
        pass

class lamina:
    def __init__(self,fiber: Material, matrix: Material, V_f: float) -> None:
        self.E_1 = calculate_E1(fiber.E,matrix.E,V_f)
        self.E_2 = calculate_E2(fiber.E,matrix.E,V_f)
        self.nu_12 = calculate_nu12(fiber.nu,matrix.nu,V_f)
        self.G_12 = calculate_G_12(fiber.G,matrix.G,V_f)
        self.matrixQ = assemble_matrixQ(self.E_1,self.E_2,self.nu_12,self.G_12) 
        pass