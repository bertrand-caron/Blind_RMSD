from yaml import dump

on_indexes = lambda chemical_point:chemical_point.index
on_elements = lambda chemical_point:chemical_point.element
on_coords = lambda chemical_point: chemical_point.x
on_canonical_rep = lambda chemical_point: chemical_point.canonical_rep


class ChemicalPoint:
    
    def __init__(self, x, index, element=None, flavour=None):
        self.x = x
        self.index = index
        self.element = element
        self.flavour = flavour
        self.canonical_rep = str(element) if self.flavour else str(element)
    
    def __str__(self):
        return '{data}'.format(data=dump({'x':self.x,
                                        'index':self.index,
                                        'element': self.element,
                                        'flavour': self.flavour,
                                        'canonical_rep': self.canonical_rep
                                        }))
    def __repr__(self):
        return self.__str__()