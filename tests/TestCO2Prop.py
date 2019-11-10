#! /bin/env python2

from __future__ import division
import unittest

from solartherm import simulation
from solartherm import postproc

from math import pi

class TestCO2Prop(unittest.TestCase):
	def setUp(self):
		fn = 'TestCO2Prop.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model()
		sim.compile_sim(args=['-s'])
		sim.simulate(start=0, stop=2, step=0.1)
		self.res = postproc.SimResult(sim.model + '_res.mat')

	def test_CO2_props(self):
		delta = 1 # error tolerated (percentage)

		self.assertLessEqual(self.res.interpolate('err_avg_h_s_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_h_T_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_cp_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_cv_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_eta_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_lambda_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_rho_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_s_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_T_h_test', 1), delta)
		self.assertLessEqual(self.res.interpolate('err_avg_U_h_test', 1), delta)
		
        

if __name__ == '__main__':
	unittest.main()
