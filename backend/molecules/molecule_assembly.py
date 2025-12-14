from rdkit import Chem

class Atom(): 

    def __init__(self, id, element, atomic_number, mass,
                 formal_charge=0, velocity=[0,0,0], position=(0,0,0)):
        
        self.id = id
        self.element = Chem.AtomFromSmarts(element)
        self.atomic_number = atomic_number
        self.valence = self.element.GetValence(Chem.ValenceType.EXPLICIT)
        self.mass = mass

        self.velocity = velocity
        self.position = position
        self.formal_charge = formal_charge

        self.hybridization = (self.valence // 2)
    
    def update_hybridization(self, bonds=(1, 1, 1, 1)):

        num_bonds = len(bonds)
        self.hybridization = (num_bonds) + ((self.valence - num_bonds) // 2)

class Bond(): 

    def __init__(self, atom_1, atom_2, order):

        self.atom_1 = atom_1
        self.atom_2 = atom_2

        self.order = order
        self.equilibrium_length = Engine.compute() 
        self.k_force_constant = Engine.compute()

        self.rotatable = True
        if order > 1: 
            self.rotatable = False

    def update_rotatable(self, is_rotatable, new_order=None):

        if new_order == 1:
            self.rotatable = is_rotatable
        else: 
            self.rotatable = is_rotatable if (self.order == 1 and new_order is None) else False

class Molecule(object): 

    def __init__(self, atoms, bonds): 

        self.atoms = atoms
        self.bonds = bonds

        # construct adjacency list 

        self.masses = []
        for atom in atoms: 
            self.masses.append(atom.mass)
        self.mol_weight = sum(self.masses)

        self.velocities = []
        for atom in atoms: 
            self.velocities.append(atom.velocity)
        self.positions = []
        for atom in atoms: 
            self.positions.append(atom.position)

        self.angles = {}
        for atom in atoms: 
            self.update_angles()
    
    def add_atom(self, id, element, atomic_number, mass, 
                formal_charge=0, velocity=[0,0,0], position=(0,0,0)):
        
        new_atom = Atom(id, element, atomic_number, mass, formal_charge, 
                            velocity, position)
        self.atoms.append(new_atom)

    def add_bond(self, atom_1, atom_2, bond_order): 

        new_bond = Bond(atom_1, atom_2, bond_order)
        self.bonds.append(new_bond)

        self.update_angles(atom_1)
        self.update_angles(atom_2)

    def update_angles(self, atom_pivot): 

        return