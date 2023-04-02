#macromechanics module for myComposite Package create for ME227 class

import numpy as np

class LaminateLayupError(TypeError): pass

def assemble_matrixT(angle):
    """ Assembles T matrix (angles transformation matrix LCS -> MCS). """

    if not isinstance(angle, (int, float)) or not (-360 <= angle <= 360):
        raise LaminateLayupError("lamina angle is not between +- 360 degrees"+
                                " or it is not an int/float")

    #Transforms angle (degrees) to angle (radians)
    angle = np.pi*angle/180

    cos = np.cos(angle)
    sin = np.sin(angle)
    cs = cos*sin
    cc = cos**2
    ss = sin**2

    T = np.zeros((3, 3))
    T = np.array([[cc,    ss,   2*cs   ],
                     [ss,    cc,  -2*cs   ],
                     [-cs, cs, cc-ss]])
    return T


# def assemble_matrixQ (mat_prop, fail_type = None):
#     """ Assembles Q matrix (reduced elastic matrix) for a given layer. """
    
#     if not isinstance(mat_prop, dict):
#         raise TypeError('mat_prop must be a dictionary')

#     # Degradation Factor (for failed layers)
#     df = 0.001

#     if fail_type == "fiber" or fail_type == "shear":
#         E1  = mat_prop["E1"]*df
#         E2  = mat_prop["E2"]*df
#         n12 = mat_prop["n12"]*df
#         G12 = mat_prop["G12"]*df
#         n21 = n12*E2/E1
#     elif fail_type == "matrix":
#         E1  = mat_prop["E1"]
#         E2  = mat_prop["E2"]*df
#         n12 = mat_prop["n12"]*df
#         G12 = mat_prop["G12"]*df
#         n21 = n12*E2/E1
#     else:
#         E1  = mat_prop["E1"]
#         E2  = mat_prop["E2"]
#         n12 = mat_prop["n12"]
#         G12 = mat_prop["G12"]
#         n21 = n12*E2/E1
    
#     Q11 = E1/(1 - n12*n21)        
#     Q12 = n12*E1*E2 / (E1 - (n12 ** 2) * E2)
#     Q22 = E1*E2 / (E1 - (n12 ** 2) * E2)
#     Q66 = G12

#     Q = np.zeros((3, 3))
#     Q = np.array([[Q11, Q12, 0],
#                      [Q12, Q22, 0],
#                      [0,   0, Q66]])
#     return Q