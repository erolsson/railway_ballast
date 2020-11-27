class BoundaryCondition(object):
    def __init__(self, set_name, bc_type, component, values=None):
        self.set_name = set_name
        self.type = bc_type
        self.component = component
        self.values = values

