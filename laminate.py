#Laminate module for myComposite Package create for ME227 class
#This module cover macromechanics calculation of a laminate
#based on classical laminate theory

from lamina import *

@dataclass
class Laminate_desc:
    """Class for containing description of a laminate. To override ply thickness value, specify thickness tuple."""
    plies: tuple[Lamina, ...]
    angles: tuple[float, ...]
    thickness: tuple[float, ...] = field(default=None)
    
    def __post_init__(self):
        if self.thickness is None:
            self.thickness = tuple(ply.thickness for ply in self.plies)

    def is_valid(self) -> bool: 
        _validity = False
        if len(self.plies) == len(self.thickness) and len(self.plies) == len(self.angles):
            _validity = True
        return _validity
    

def assemble_Z_vector(desc: Laminate_desc):
    """ Assembles Z vector which contains z coordinates of laminate (with notation described in  Jone's Book). """
    if not desc.is_valid():
        print ("catched expected error")
        raise Exception("The Laminate description is not valid")
    else:
        layer_count = len(desc.thickness)
        _total_thickness = np.sum(desc.thickness)
        _Z = np.zeros(layer_count+1, dtype=float)
        _Z[0] = -_total_thickness/2.0

        for i in range(layer_count):
            _Z[i+1] = _Z[i] + desc.thickness[i]

    return _Z

def assemble_ABD_matrix(desc: Laminate_desc):
    """ Assembles ABD matrix (laminate stiffness matrix). """
    if not desc.is_valid():
        print ("catched expected error")
        raise Exception("The Laminate description is not valid")

    _Z = assemble_Z_vector(desc)
    layer_count = len(desc.thickness)
    _A = np.zeros((3,3))
    _B = np.zeros((3,3))
    _D = np.zeros((3,3))
    for k in range(layer_count):
        _Qbar = calculate_matrixQbar(desc.plies[k].matrixQ, desc.angles[k])
        #assemble A
        _A = _A + _Qbar*(_Z[k+1]-_Z[k])
        #assemble B
        _B = _B + (1/2)*_Qbar*(_Z[k+1]**2 - _Z[k]**2) #This need - to match result from online solver.
        #assemble D
        _D = _D + (1/3)*_Qbar*(_Z[k+1]**3 - _Z[k]**3)
    
    _ABD = np.block([
        [_A,_B],
        [_B,_D]
    ])

    return _ABD

class Laminate:
    def __init__(self,desc: Laminate_desc,name:str = "unnamed") -> None:
        self.name = name
        self.description = desc
        self.vectorZ = assemble_Z_vector(self.description)
        self.matrixABD = assemble_ABD_matrix(self.description)
        # self.global_stiffness_matrix = np.sum([ply.matrixQ* for ply in desc.plies],axis=0)
        self.density = np.sum([ply.density*ply.thickness for ply in desc.plies])/np.sum([ply.thickness for ply in desc.plies])
        pass

    def get_matrixA(self):
        return self.matrixABD[:3,:3]
    
    def get_matrixB(self):
        return self.matrixABD[:3,-3:]
    
    def get_matrixD(self):
        return self.matrixABD[-3:,-3:]