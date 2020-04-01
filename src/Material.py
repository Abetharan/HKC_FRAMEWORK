class Material:

    def __init__(self, z, atomicnumber, idealgas = True, feospath1 = None):
        
        self.z = z
        self.ar = atomicnumber
        self.idealgas = idealgas
        self.feospath1 = feospath1