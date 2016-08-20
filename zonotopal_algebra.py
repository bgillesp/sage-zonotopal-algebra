class ZonotopalAlgebra(object):
    def __init__(self,X,variant="central",data={}):
        self.variant = variant
        self.data = data
        self.obj = None
        if type=="central":
            varNames = 'x'; # TODO default variable
            if 'varNames' in data:
                varNames = data['varNames']
            self.obj = new CentralZonotopalAlgebra(X, varNames)
        # TODO include additional types of zonotopal algebras
        else:
            raise ValueException("uncrecognized zonotopal algebra type: %s" % type)

    def I_ideal_gens(self):
        return self.obj.I_ideal_gens()

    def J_ideal_gens(self):
        return self.obj.I_ideal_gens()

    def D_basis(self):
        return self.obj.D_basis()

    def P_basis(self):
        return self.obj.P_basis()


    
