import numpy as np
from scipy.integrate import solve_ivp


def invariant_1(tensor):
    return tensor[0] + tensor[1] + tensor[2]


def von_mises(tensor):
    return np.sqrt(np.sum(tensor[0:3]**2) + 3*np.sum(tensor[3:6]**2) - tensor[0]*tensor[1]
                   - tensor[0]*tensor[2] - tensor[1]*tensor[2])


class MaterialModel:
    def __init__(self, material_parameters):
        self.gf = abs(material_parameters[0])
        self.A = abs(material_parameters[1])
        self.A1 = abs(material_parameters[2])
        self.A2 = material_parameters[3]

        self.nf = abs(material_parameters[4])
        self.H1 = abs(material_parameters[5])
        self.b1 = abs(material_parameters[6])
        self.b2 = np.exp(-abs(material_parameters[7]))
        self.b3 = abs(material_parameters[8])
        self.nb = abs(material_parameters[9])
        self.b4 = abs(material_parameters[10])
        self.b5 = abs(material_parameters[11])
        self.b6 = material_parameters[12]

        self.frictional_strain = None
        self.compaction_strain = None
        self.parameters = material_parameters

    def update(self, cycles, cyclic_stress, static_stress):
        self.frictional_strain = np.zeros((cycles.shape[0], 6))
        self.compaction_strain = np.zeros((cycles.shape[0], 6))
        p0 = -invariant_1(static_stress)/3
        pc = -invariant_1(cyclic_stress)/3
        if p0 < 0:
            p0 = 0

        q = von_mises(cyclic_stress)
        nij = 1.5*(cyclic_stress - invariant_1(cyclic_stress)/3*np.array([1, 1, 1, 0, 0, 0]))/q

        def dkd_dn(_, ep):
            arg = 1. + self.A1*p0 + self.A2*p0**2
            hf = self._hf(von_mises(ep))

            if arg < 1e-6:
                arg = 1e-6
            g = np.sqrt(arg)
            f = (q/g - hf)
            if f < 0.:
                f = 0.

            ep_eff_dn = self.A*f**self.gf
            dilatation = self.b1 - self.b5*p0 - self.b6*p0**2
            deij_dn = ep_eff_dn*(nij + dilatation*np.array([1, 1, 1, 0, 0, 0]))
            return deij_dn

        def dkv_dn(_, ev):

            f = (self.b4*p0 + pc - self.b3*ev)
            if f < 0:
                f = 0
            return self.b2*f**self.nb

        for i in range(1, cycles.shape[0]):
            e_fric = self.frictional_strain[i - 1, :]
            e_comp = self.compaction_strain[i - 1, :]
            n0 = cycles[i - 1]
            solution_fric = solve_ivp(dkd_dn, [n0, cycles[i]], e_fric)
            solution_comp = solve_ivp(dkv_dn, [n0, cycles[i]], [e_comp[0] + e_comp[1] + e_comp[2]])

            e_f = solution_fric.y[:, -1]
            e_c = solution_comp.y[:, -1]
            if np.linalg.norm(e_f) > 1.:
                e_f /= np.linalg.norm(e_f)
            self.frictional_strain[i, :] = e_f
            for j in range(3):
                if abs(e_c) > 1:
                    e_c /= abs(e_c)
                self.compaction_strain[i, j] = e_c/3

        return self.frictional_strain - self.compaction_strain

    def _hf(self, ep):
        return self.H1*(1 - np.exp(-self.nf*ep))

    def strain(self):
        e = self.frictional_strain - self.compaction_strain
        e[:, 3:] *= 2
        return e

    def volumetric_strain(self):
        e = self.strain()
        return e[:, 0] + e[:, 1] + e[:, 2]

    def deviatoric_strain(self):
        e_dev = np.array(self.strain())
        e_dev[:, 3:] *= 2
        e_vol = self.volumetric_strain()
        for i in range(3):
            e_dev[:, i] -= e_vol/3
        return e_dev
