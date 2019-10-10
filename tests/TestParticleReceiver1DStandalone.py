#! /bin/env python2
from __future__ import division
import unittest
from solartherm import simulation, postproc

class TestParticleReceiver1DStandalone(unittest.TestCase):
	def test1(self):
		print "RUNNING SETUP"
		fn = 'ParticleReceiver1D_standalone.mo'
		sim = simulation.Simulator(fn)
		print "COMPILING MODEL"
		sim.compile_model()
		sim.compile_sim(args=['-s'])
		print "SOLVING MODEL"
		sim.simulate(start=0, stop='1s', step='1s', solver='dassl',nls=None)
		self.res = postproc.SimResult(sim.model + '_res.mat')

		print "RESULTS"
		N = self.res.interpolate('N',1)
		vl = ['N','T_s[1]','T_s[%d]'%(N,), 'm_flow', 'T_a']
		for v in vl:
			print '%s = %f' %(v,self.res.interpolate(v,1))

if __name__ == '__main__':
	unittest.main()
