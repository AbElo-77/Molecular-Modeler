import numpy as np 

class Angle(): 

    def __init__(self, atom_1, atom_2, atom_3):

        self.vector_1_2 = (atom_2.position) - (atom_1.position)
        self.vector_2_3 = (atom_3.position) - (atom_2.position)

        cosine_theta = np.dot(self.vector_1_2, self.vector_1_2)
        cosine_theta = cosine_theta / (np.linalg.norm(self.vector_1_2) * np.linalg.norm(self.vector_1_2))
        self.theta = np.arccos(cosine_theta)

