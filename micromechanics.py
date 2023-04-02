#micromechanics module for myComposite Package create for ME227 class
#This module cover calculation of a lamina properties using 

import numpy as np

class fiber:
    def __init__(self,name: str, modulus: float, poisson_ratio:float) -> None:
        self.name = name
        self.modulus = modulus
        self.poisson_ratio = poisson_ratio
        pass

