import numpy as np
from scipy.integrate import solve_ivp

from experimental_results import sun_et_al_16


def invariant_1(tensor):
    return tensor[0] + tensor[1] + tensor[2]


def von_mises(tensor):
    return np.sqrt(np.sum(tensor[0:3]**2) + 3*np.sum(tensor[3:6]**2) - tensor[0]*tensor[1]
                   - tensor[0]*tensor[2] - tensor[1]*tensor[2])


class MaterialModel:
    def __init__(self, parameters, frequency=5.):
        self.gf = abs(parameters[0])
        self.A = abs(parameters[1])
        self.A1 = abs(parameters[2])
        self.A2 = parameters[3]

        self.nf = abs(parameters[4])
        self.H1 = abs(parameters[5])
        self.b1 = parameters[9]
        self.b2 = parameters[10]
        self.b3 = parameters[11]
        self.nb = parameters[12]

        freq_idx = {10.: 6, 20.: 7, 40.: 8}
        if frequency == 5.:
            self.fd = 1.
        else:
            self.fd = parameters[freq_idx[frequency]]
            if self.fd < 1.:
                self.fd = 1
        self.strain = None
        self.parameters = parameters

    def update(self, cycles, cyclic_stress, static_stress):
        self.strain = np.zeros((cycles.shape[0], 6))
        p = -invariant_1(static_stress)/3
        if p < 0:
            p = 0
        q = von_mises(cyclic_stress)
        p_cyclic = -invariant_1(cyclic_stress)/3
        p_eff = p + p_cyclic
        nij = 1.5*(cyclic_stress - invariant_1(cyclic_stress)/3*np.array([1, 1, 1, 0, 0, 0]))/q

        def dedn(_, ep):
            arg = 1. + self.A1*p_eff + self.A2*p_eff**2
            hf = self._hf(von_mises(ep))

            if arg < 1e-6:
                arg = 1e-6
            g = np.sqrt(arg)
            f = (self.fd*q/g - hf)
            if f < 0.:
                f = 0.

            ep_eff_dn = self.A*f**self.gf
            dilatation = self.b1*np.exp(-self.nf*von_mises(ep)) + self.b2*p + self.b3*p_cyclic**self.nb
            deij_dn = ep_eff_dn*(nij + dilatation*np.array([1, 1, 1, 0, 0, 0]))
            # e = ep_eff_dN*(nij + (self.b1 + self.b2*(p + self.fd*p_cyclic))*np.array([1., 1., 1., 0., 0., 0]))

            return deij_dn

        for i in range(cycles.shape[0]):
            if i == 0:
                n0 = 0
                e1 = np.zeros(6)
            else:
                e1 = self.strain[i - 1, :]
                n0 = cycles[i - 1]

            solution = solve_ivp(dedn, [n0, cycles[i]], e1)

            for j in range(6):
                e = solution.y[j][-1]
                if e > 1.:
                    e = 1
                if e < -1:
                    e = -1
                self.strain[i, j] = e
        return self.strain

    def _hf(self, ep):
        return self.H1*(1 - np.exp(-self.nf*ep))

    def permanent_strain_tensor(self):
        return self.strain

    def volumetric_strain(self):
        eij = self.permanent_strain_tensor()
        ev = np.sum(eij[:, 0:3], axis=1)
        # ev[np.abs(ev) > 0.25] = -1.
        return ev

    def deviatoric_strain(self):
        eij = self.permanent_strain_tensor()
        e_vol = self.volumetric_strain()
        edev = eij - np.outer(e_vol/3, np.array([1, 1, 1, 0, 0, 0]))
        # edev[edev > 0.35] = 1.
        # edev[edev < -0.35] = -1.
        return edev


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
    par = np.array([1.92596342e+00, 1.33771787e-07, 2.36556786e-02, 6.54713185e-04,
                    1.59527656e+01, 2.03974910e+02, 1., 1.,
                    1., 9.05538667e-02, -6.54359191e-03,  7.15099017e-07,  2.62519248e+00])

    cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))

    for i, f in enumerate([5., 10., 20., 40.]):
        model = MaterialModel(parameters=par, frequency=f)
        experimental_data = sun_et_al_16.get_data(f=f)
        for experiment in experimental_data:
            plt.figure(i)
            edev = experiment.deviatoric_axial_strain()
            plt.semilogx(experiment.cycles[edev < 0.9], edev[edev < 0.9])
            static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
            cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
            model.update(cycles, cyclic_stress, static_stress)
            edev = -model.deviatoric_strain()[:, 0]
            plt.semilogx(cycles[edev < 0.9], edev[edev < 0.9] + experiment.deviatoric_axial_strain()[0], '--')

            plt.figure(i+4)
            plt.semilogx(experiment.cycles, experiment.volumetric_strain)
            plt.semilogx(cycles, -model.volumetric_strain(), '--')

        plt.figure(i)
        plt.xlim(1, 5e5)
        plt.figure(i+4)
        plt.xlim(1, 5e5)

    plt.figure(100)
    p = np.linspace(-100, 100, 100)
    plt.plot(p, np.sqrt(1 + par[2]*p + par[3]*p**2))
    plt.show()


if __name__ == '__main__':
    main()
