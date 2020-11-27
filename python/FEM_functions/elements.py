from __future__ import print_function

from collections import namedtuple

# import find_modules  # noqa
import numpy as np


class C3D8:
    local_nodal_pos = np.array([[-1, -1, -1],
                                [1, -1, -1],
                                [1, 1, -1],
                                [-1, 1, -1],
                                [-1, -1, 1],
                                [1, -1, 1],
                                [1, 1, 1],
                                [-1, 1, 1]])
    gauss_points = np.zeros((8, 3))
    counter = 0
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                gauss_points[counter, :] = k/np.sqrt(3), j/np.sqrt(3), i/np.sqrt(3)
                counter += 1

    def __init__(self, nodes):
        self.xe = np.zeros((8, 3))
        for i in range(3):
            self.xe[:, i] = [n.coordinates[i] for n in nodes]
            self.node_labels = [n.label for n in nodes]

    def d(self, xi, eta, zeta):
        """
        Returns the a matrix 3x8 matrix with derivatives of the shape functions with respect to x y and z
        :param xi:     xi-coordinate
        :param eta:    eta-coordinate
        :param zeta:   zeta-coordinate
        :return:       3x8 matrix with derivatives
        """

        d_matrix = np.zeros((3, 8))
        for i in range(8):
            d_matrix[0, i] = (1. + eta*self.local_nodal_pos[i, 1]) * \
                             (1. + zeta*self.local_nodal_pos[i, 2])*self.local_nodal_pos[i, 0]/8
            d_matrix[1, i] = (1. + xi*self.local_nodal_pos[i, 0]) * \
                             (1. + zeta*self.local_nodal_pos[i, 2])*self.local_nodal_pos[i, 1]/8
            d_matrix[2, i] = (1. + xi*self.local_nodal_pos[i, 0]) * \
                             (1. + eta*self.local_nodal_pos[i, 1])*self.local_nodal_pos[i, 2]/8
        return d_matrix

    def J(self, xi, eta, zeta):
        return np.dot(self.d(xi, eta, zeta), self.xe)

    def B(self, xi, eta, zeta):
        B = np.zeros((6, 24))
        jacobian = self.J(xi, eta, zeta)
        d = self.d(xi, eta, zeta)
        for i in range(8):
            dx_avg = [np.linalg.solve(self.J(*gp), self.d(*gp)[:, i])*np.linalg.det(self.J(*gp))
                      for gp in self.gauss_points]
            dx_avg = sum(dx_avg)/self.volume()

            for j in range(3):
                for k in range(3):
                    B[j, 3*i + k] += dx_avg[k]/3

            dx = np.linalg.solve(jacobian, d[:, i])
            for j in range(3):
                for k in range(3):
                    if j == k:
                        B[j, 3*i + k] += 2*dx[k]/3
                    else:
                        B[j, 3*i + k] += -dx[k]/3

            B[3, 3*i] += dx[1]
            B[3, 3*i + 1] += dx[0]

            B[4, 3*i] += dx[2]
            B[4, 3*i + 2] += dx[0]

            B[5, 3*i + 1] += dx[2]
            B[5, 3*i + 2] += dx[1]
        return B

    def volume(self):
        return np.sum([np.linalg.det(self.J(*gp)) for gp in self.gauss_points])


if __name__ == '__main__':
    Node = namedtuple('Node', ['coordinates', 'label'])
    node_list = [Node(coordinates=[0, 0, 0], label=1),
                 Node(coordinates=[2, 0, 0], label=2),
                 Node(coordinates=[2, 2, 0], label=3),
                 Node(coordinates=[0, 2, 0], label=4),
                 Node(coordinates=[0, 0, 1], label=5),
                 Node(coordinates=[2, 0, 1], label=6),
                 Node(coordinates=[2, 2, 1], label=7),
                 Node(coordinates=[0, 2, 1], label=8)]
    element = C3D8(node_list)
    element.B(0, 0, 0)
    print(np.linalg.solve(element.J(0, 0, 0), element.d(0, 0, 0)[:, 0])*8/element.volume())
