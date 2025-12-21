from rdkit import Chem
import numpy as np
from backend.molecules.geometry_logic import Angle, Dihedral

class Atom(): 

    def __init__(self, id, element,
                 formal_charge=0, velocity=None, position=None):
        
        self.id = id
        self.element = Chem.Atom(element)
        self.atomic_number = self.element.GetAtomicNum()
        self.valence = self.element.GetValence(Chem.ValenceType.EXPLICIT)
        self.mass = self.element.GetMass()

        zero_array = np.zeros(3)
        self.velocity = velocity if velocity is not None else zero_array
        self.position = position if position is not None else zero_array

        self.formal_charge = formal_charge
        self.hybridization = None
    
    def update_hybridization(self, hybridization):

        self.hybridization = hybridization

class Bond(): 

    def __init__(self, atom_1, atom_2, order):

        self.atom_1 = atom_1
        self.atom_2 = atom_2

        self.order = order
        # self.equilibrium_length = Engine.compute() 
        # self.k_force_constant = Engine.compute()

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

        self.adjacency_list = {}
        for atom in atoms: 
            adjacent_atoms = []
            for bond in bonds: 
                if bond.atom_1 == atom or bond.atom_2 == atom: 
                    adjacent_atoms.append(bond.atom_2 if bond.atom_1 == atom else bond.atom_1)
            self.adjacency_list[atom] = adjacent_atoms

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

        self.bond_angles = {}
        self.torsion_angles = {}

        for atom in atoms: 
            self.update_angles(atom)

    
    def add_atom(self, id, element,
                formal_charge=0, velocity=None, position=None):
        
        new_atom = Atom(id, element, formal_charge, 
                            velocity, position)
        self.atoms.append(new_atom)
        self.adjacency_list[new_atom] = []

    def add_bond(self, atom_1, atom_2, bond_order): 

        new_bond = Bond(atom_1, atom_2, bond_order)
        self.bonds.append(new_bond)

        self.adjacency_list[atom_1].append(atom_2)
        self.adjacency_list[atom_2].append(atom_1)

        self.update_angles(atom_1)
        self.update_angles(atom_2)

    def update_angles(self, atom_pivot): 

        adjacent_atoms = self.adjacency_list[atom_pivot]

        for atom in adjacent_atoms:
            start_index = 1
            for i in range(start_index, len(adjacent_atoms)):
                if (atom, atom_pivot, adjacent_atoms[i]) in self.bond_angles or (adjacent_atoms[i], atom_pivot, atom) in self.bond_angles:
                    continue

                params = (0, 0)
                # params = Engine.get_params(atom, atom_pivot, adjacent_atoms[i])

                angle_new = Angle(atom, atom_pivot, adjacent_atoms[i], params)
                self.bond_angles[(atom, atom_pivot, adjacent_atoms[i])] = angle_new

        for atom in adjacent_atoms:
            self.update_dihedrals_for_bond(atom_pivot, atom)
            
    def update_dihedrals_for_bond(self, atom_j, atom_k):

        neighbors_j = self.adjacency_list[atom_j]
        neighbors_k = self.adjacency_list[atom_k]

        for atom_i in neighbors_j:
            if atom_i == atom_k:
                continue

            for atom_l in neighbors_k:
                if atom_l == atom_j:
                    continue

                key = (atom_i, atom_j, atom_k, atom_l)
                rev_key = (atom_l, atom_k, atom_j, atom_i)

                if key not in self.torsion_angles and rev_key not in self.torsion_angles:

                    params = (0, 0)
                    # params = Engine.get_params(atom_i, atom_j, atom_k, atom_l)

                    self.torsion_angles[key] = Dihedral(atom_i, atom_j, atom_k, atom_l, params)