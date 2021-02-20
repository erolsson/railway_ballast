import numpy as np
from scipy.integrate import solve_ivp

from experimental_results import sun_et_al_16
from model_parameters import parameters, parameters_common


def invariant_1(tensor):
    return tensor[0] + tensor[1] + tensor[2]


def von_mises(tensor):
    return np.sqrt(np.sum(tensor[0:3]**2) + 3*np.sum(tensor[3:6]**2) - tensor[0]*tensor[1]
                   - tensor[0]*tensor[2] - tensor[1]*tensor[2])


class MaterialModel:
    def __init__(self, material_parameters, frequency=5.):
        self.gf = abs(material_parameters[0])
        self.A = abs(material_parameters[1])
        self.A1 = abs(material_parameters[2])
        self.A2 = material_parameters[3]

        self.nf = abs(material_parameters[4])
        self.H1 = abs(material_parameters[5])
        self.b1 = abs(material_parameters[9])
        self.b2 = np.exp(-abs(material_parameters[10]))
        self.b3 = abs(material_parameters[11])
        self.nb = abs(material_parameters[12])
        self.b4 = abs(material_parameters[13])
        self.b5 = abs(material_parameters[17])
        self.b6 = material_parameters[18]

        freq_idx = {10.: 6, 20.: 7, 40.: 8}
        c_idx = {10.: 14, 20.: 15, 40.: 16}
        if frequency == 5.:
            self.fd = 1.
            self.fc = 0.
        else:
            self.fd = material_parameters[freq_idx[frequency]]
            self.fc = abs(material_parameters[c_idx[frequency]])
            if self.fd > 1.:
                self.fd = 1

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
        pm = (p0 + pc)/2
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
            dilatation = self.b1 - self.b5*p0 - self.b6*p0**2 - self.fc
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

            self.frictional_strain[i, :] = e_f
            for j in range(3):
                self.compaction_strain[i, j] = e_c/3
        return self.frictional_strain - self.compaction_strain

    def _hf(self, ep):
        return self.H1*(1 - np.exp(-self.nf*self.fd*ep))

    def strain(self):
        return self.frictional_strain - self.compaction_strain

    def volumetric_strain(self):
        e = self.strain()
        return e[:, 0] + e[:, 1] + e[:, 2]

    def deviatoric_strain(self):
        e_dev = np.array(self.strain())
        e_vol = self.volumetric_strain()
        for i in range(3):
            e_dev[:, i] -= e_vol/3
        return e_dev


def main():
    import matplotlib.pyplot as plt
    import matplotlib.style

    matplotlib.style.use('classic')
    plt.rc('text', usetex=True)
    plt.rc('font', serif='Computer Modern Roman')
    plt.rcParams.update({'font.size': 20})
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                      'monospace': ['Computer Modern Typewriter']})
    # 5 Hz
    base_parameters = np.zeros(19)
    base_parameters[6:9] = 1
    base_parameters[14:17] = 0
    # base_parameters[:9] = parameters_common
    cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))

    colors = {(10, 230): 'b', (30, 230): 'r', (60, 230): 'g', (60, 370): 'k', (30, 276): 'm', (60, 460): 'y'}
    for i, f in enumerate([5., 10., 20., 40.]):
        # par[0:6] = parameters[f]

        experimental_data = sun_et_al_16.get_data(f=f)
        print("\nf=" + str(f) + " Hz")
        par1 = np.array(base_parameters)
        par1[0:6] = parameters[f][0:6]
        par1[9:14] = parameters[f][6:11]
        par1[17:19] = parameters[f][11:13]
        print(par1)
        for experiment in experimental_data:
            p = experiment.p
            q = experiment.q
            plt.figure(i)
            edev = experiment.deviatoric_axial_strain()
            plt.semilogx(experiment.cycles[edev < 0.9], edev[edev < 0.9], colors[(p, q)], lw=2)

            model = MaterialModel(material_parameters=par1, frequency=f)
            static_stress = -p*np.array([1, 1, 1, 0, 0, 0])
            cyclic_stress = -q*np.array([1, 0, 0, 0, 0, 0])
            model.update(cycles, cyclic_stress, static_stress)
            edev = -model.deviatoric_strain()[:, 0]
            plt.semilogx(cycles[edev < 0.9], edev[edev < 0.9] + experiment.deviatoric_axial_strain()[0],
                         '--' + colors[(p, q)], lw=2)

            plt.figure(i + 4)
            plt.semilogx(experiment.cycles, experiment.volumetric_strain, colors[(p, q)], lw=2)
            plt.semilogx(cycles, -model.volumetric_strain() + experiment.volumetric_strain[0],
                         '--' + colors[(p, q)], lw=2)

            plt.figure(i)
            par2 = parameters_common
            model = MaterialModel(material_parameters=par2, frequency=f)
            model.update(cycles, cyclic_stress, static_stress)
            edev = -model.deviatoric_strain()[:, 0]
            plt.semilogx(cycles[edev < 0.9], edev[edev < 0.9] + experiment.deviatoric_axial_strain()[0],
                         ':' + colors[(p, q)], lw=2)
            plt.figure(i + 4)
            plt.semilogx(cycles, -model.volumetric_strain() + experiment.volumetric_strain[0],
                         ':' + colors[(p, q)], lw=2)
            print(experiment.volumetric_strain[-1],
                  -model.volumetric_strain()[-1] + experiment.volumetric_strain[0],
                  experiment.volumetric_strain[0])
        plt.figure(i)

        plt.xlim(1, 5e5)

    plt.show()


if __name__ == '__main__':
    main()
