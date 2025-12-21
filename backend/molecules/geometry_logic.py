import numpy as np

class Angle:

    def __init__(self, atom_1, atom_2, atom_3, params):

        self.theta0, self.k_constant = params

        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self.atom_3 = atom_3

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
            return []

        normal_vector /= norm
        return normal_vector
    
    def energy(self): 

        theta = self.calculate_value()
        return 0.5 * self.k_constant * (theta - self.theta0)**2
    
    def forces(self):

        dot = np.dot(self.vector_1_2, self.vector_2_3)
        norm = (np.linalg.norm(self.vector_1_2) * np.linalg.norm(self.vector_2_3))

        cosine_theta = dot / norm
        cosine_theta = np.clip(cosine_theta, -1.0, 1.0)

        dU = self.k_constant * (self.calculate_value() - self.theta0)
        sin_theta = np.sqrt(1 - cosine_theta**2)

        if sin_theta < 1e-12: 
            return {}
        force_coefficent = -dU / sin_theta

        force_on_1 = force_coefficent * (self.vector_2_3 / np.linalg.norm(self.vector_2_3) 
                                         - cosine_theta *  self.vector_1_2  / np.linalg.norm(self.vector_1_2)) / np.linalg.norm(self.vector_1_2)
        force_on_2 = force_coefficent * (self.vector_1_2 / np.linalg.norm(self.vector_1_2) 
                                         - cosine_theta * self.vector_2_3  / np.linalg.norm(self.vector_2_3)) / np.linalg.norm(self.vector_2_3)
        force_on_3 = -(force_on_1 + force_on_2)      

        return {
            self.atom_1: force_on_1, 
            self.atom_2: force_on_3, 
            self.atom_3: force_on_2
        }

class Dihedral:

    def __init__(self, atom_1, atom_2, atom_3, atom_4, params):

        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self.atom_3 = atom_3
        self.atom_4 = atom_4

        self.b1 = atom_2.position - atom_1.position
        self.b2 = atom_3.position - atom_2.position
        self.b3 = atom_4.position - atom_3.position
        self.params = params

        self.plane_1_2_3 = np.cross(self.b1, self.b2)
        self.plane_2_3_4 = np.cross(self.b2, self.b3)
        
    def calculate_value(self): 

        n1_norm = np.linalg.norm(self.plane_1_2_3)
        n2_norm = np.linalg.norm(self.plane_2_3_4)

        self.plane_1_2_3 /= n1_norm
        self.plane_2_3_4 /= n2_norm

        b2_hat = self.b2 / np.linalg.norm(self.b2)
        x = np.dot(self.plane_1_2_3, self.plane_2_3_4)
        y = np.dot(b2_hat, np.cross(self.plane_1_2_3, self.plane_2_3_4))

        self.phi = np.arctan2(y, x)

    def energy(self):

        phi = self.calculate_value()
        U = 0.0
        for n, gamma, Vn in self.params:
            U += 0.5 * Vn * (1 + np.cos(n*phi - gamma))
        return U
    
    def forces(self):
        
        phi = self.calculate_value()
        torque = 0
        for n, gamma, Vn in self.params:
            torque += -0.5 * Vn * n * np.sin(n*phi - gamma)

        force_on_1 = -torque * (self.plane_1_2_3 / np.linalg.norm(self.plane_1_2_3)) / np.linalg.norm(self.b1)**2
        force_on_4 = -torque * (self.plane_2_3_4 / np.linalg.norm(self.plane_2_3_4)) / np.linalg.norm(self.b3)**2

        alpha = np.dot(force_on_1 + force_on_4, self.b2) / np.linalg.norm(self.b2)**2

        force_on_2 = -force_on_1 + alpha * self.b2
        force_on_3 = -force_on_4 - alpha * self.b2

        total_force = force_on_1 + force_on_2 + force_on_3 + force_on_4
        assert np.allclose(total_force, 0, atol=1e-10), "Conservation Law Violated! Oh No"

        return {
            self.atom_1: force_on_1,
            self.atom_2: force_on_2,
            self.atom_3: force_on_3,
            self.atom_4: force_on_4
        }

