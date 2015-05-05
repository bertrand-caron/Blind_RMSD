from yaml import dump

on_indexes = lambda chemical_point:chemical_point.index
on_elements = lambda chemical_point:chemical_point.element
on_coords = lambda chemical_point: chemical_point.x
on_canonical_rep = lambda chemical_point: chemical_point.canonical_rep


class ChemicalPoint:
    
    def __init__(self, x, index, element=None, flavour=None, grouped_flavours=None):
        self.x = x
        self.index = index
        self.element = element
        self.flavour = flavour
        self.canonical_rep = element if not grouped_flavours else '{element}{group_length}'.format(element=element, group_length=len(grouped_flavours[flavour]))
    
    def __str__(self):
        return '{{index={index}, element={element}, canonical_rep={canonical_rep}}}'.format(
                                        **{   'index':self.index,
                                            'element': self.element,
                                            'canonical_rep': self.canonical_rep
                                        })
    def __repr__(self):
        return self.__str__()
