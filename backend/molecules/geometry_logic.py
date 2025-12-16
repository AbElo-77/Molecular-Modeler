import numpy as np

class Angle:

    def __init__(self, atom_1, atom_2, atom_3, params):

        self.theta0, self.k_constant = params

        self.vector_1_2 = atom_1.position - atom_2.position
        self.vector_2_3 = atom_3.position - atom_2.position

    def calculate_value(self): 

        dot = np.dot(self.vector_1_2, self.vector_2_3)
        norm = (np.linalg.norm(self.vector_1_2) * np.linalg.norm(self.vector_2_3))

        cosine_theta = dot / norm
        cosine_theta = np.clip(cosine_theta, -1.0, 1.0)

        return np.arccos(cosine_theta)

    def get_plane(self):

        normal_vector = np.cross(self.vector_1_2, self.vector_2_3)
        norm = np.linalg.norm(normal_vector)

        if norm < 1e-12:
            return None

        normal_vector /= norm
        return normal_vector
    
    def energy(self): 

        theta = self.calculate_value()
        return 0.5 * self.k * (theta - self.theta0)**2

class Dihedral:

    def __init__(self, atom_1, atom_2, atom_3, atom_4):

        b1 = atom_2.position - atom_1.position
        b2 = atom_3.position - atom_2.position
        b3 = atom_4.position - atom_3.position

        self.plane_1_2_3 = np.cross(b1, b2)
        self.plane_2_3_4 = np.cross(b2, b3)

        n1_norm = np.linalg.norm(self.plane_1_2_3)
        n2_norm = np.linalg.norm(self.plane_2_3_4)

        self.plane_1_2_3 /= n1_norm
        self.plane_2_3_4 /= n2_norm

        b2_hat = b2 / np.linalg.norm(b2)
        x = np.dot(self.plane_1_2_3, self.plane_2_3_4)
        y = np.dot(b2_hat, np.cross(self.plane_1_2_3, self.plane_2_3_4))

        self.phi = np.arctan2(y, x)
