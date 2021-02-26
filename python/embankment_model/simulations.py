from collections import namedtuple
from copy import deepcopy

import numpy as np


ConcreteSlab = namedtuple('ConcreteSlab', ['width', 'height'])


class Sleeper:
    def __init__(self, cc_distance, width, length, height, center_sleeper=True):
        self.cc_distance = cc_distance
        self.width = width
        self.length = length
        self.height = height
        self.center_sleeper = center_sleeper


class ElasticLayer:
    def __init__(self, name, bottom_width, top_width, height, E, v, density):
        self.name = name
        self.bottom_width = bottom_width
        self.top_width = top_width
        self.height = height
        self.E = E
        self.v = v
        self.density = density


class Simulation:
    def __init__(self, job_name, cycles, layers, embankment_length, axes_load):
        self.job_name = job_name
        self.cycles = cycles
        self.layers = layers
        self.embankment_length = embankment_length
        self.track_gauge = 1.435
        self.axes_load = axes_load
        self.concrete_slab = None
        self.sleepers = None

    @property
    def track_lower_height(self):
        h = sum([layer.height for layer in self.layers])
        if self.concrete_slab:
            return h + self.concrete_slab.height
        elif self.sleepers:
            return h + self.sleepers.height
        return h


sleepers = Sleeper(cc_distance=0.65, width=0.265, length=1., height=0.175)
concrete_slab = ConcreteSlab(width=1.25, height=0.5)

simulation1 = Simulation(job_name='embankment_second_order_2', cycles=np.array([100, 1000, 10000, 100000, 1000000]),
                         embankment_length=5.5, axes_load=22.5,
                         layers=[ElasticLayer(name='soil', bottom_width=20, top_width=20.,
                                              height=2., E=50e6, v=0.2, density=2000),
                                 ElasticLayer(name='clay', bottom_width=20, top_width=20.,
                                              height=4., E=50e6, v=0.2, density=2000),
                                 ElasticLayer(name='subgrade', bottom_width=20, top_width=20.,
                                              height=1., E=50e6, v=0.2, density=2000),
                                 ElasticLayer(name='ballast_1', bottom_width=9.3, top_width=3.6,
                                              height=3.8, E=200e6, v=0.35, density=1600),
                                 ElasticLayer(name='ballast_2', bottom_width=2.55, top_width=1.8,
                                              height=0.5, E=200e6, v=0.35, density=1600)])

embankment_21_22_5t = Simulation(job_name='embankment_21_22_5t', cycles=np.array([100, 1000, 10000, 100000, 1000000]),
                                 embankment_length=5.5, axes_load=22.5,
                                 layers=[ElasticLayer(name='soil', bottom_width=20, top_width=20.,
                                                      height=2., E=50e6, v=0.2, density=2000),
                                         ElasticLayer(name='clay', bottom_width=20, top_width=20.,
                                                      height=4., E=50e6, v=0.2, density=2000),
                                         ElasticLayer(name='subgrade', bottom_width=20, top_width=20.,
                                                      height=1., E=50e6, v=0.2, density=2000),
                                         ElasticLayer(name='ballast_1', bottom_width=11.2, top_width=3.6,
                                                      height=3.8, E=200e6, v=0.35, density=1600),
                                         ElasticLayer(name='ballast_2', bottom_width=2.55, top_width=1.8,
                                                      height=0.5, E=200e6, v=0.35, density=1600)])

simulation1.sleepers = sleepers

simulation2 = deepcopy(simulation1)
embankment_21_22_5t_slab = deepcopy(embankment_21_22_5t)
embankment_21_22_5t.sleepers = sleepers
embankment_21_22_5t_slab.job_name += '_slab'
embankment_21_22_5t_slab.concrete_slab = concrete_slab
simulation2.axes_load = 30.
simulation2.job_name = 'embankment_30t'
