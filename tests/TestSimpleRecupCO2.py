#! /bin/env python2

from __future__ import division
import unittest

from solartherm import simulation
from solartherm import postproc

from math import pi

class TestSimpleRecupCO2(unittest.TestCase):
	def setUp(self):
		fn = 'TestSimpleRecupCO2.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model()
		sim.compile_sim(args=['-s'])
		sim.simulate(start=0, stop=2, step=0.1)
		self.res = postproc.SimResult(sim.model + '_res.mat')

	def test_moltenSalt_props(self):
		delta = 1 # error tolerated (percentage)

		self.assertLessEqual(self.res.interpolate('err_m_flow', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_h_HTR_turb', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_h_HTR_comp', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_deltaT_HTR', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_h_in_comp', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_Q_cooler', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_Q_heater', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_Q_HTR', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_eta', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_W_comp', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_W_turb', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_P', 1), delta)

if __name__ == '__main__':
	unittest.main()
