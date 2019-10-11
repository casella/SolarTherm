#! /bin/env python2
from __future__ import division
import unittest
from solartherm import simulation, postproc
import os.path

PLOTME=0
RUNSIM=1

class TestParticleReceiver1DStandalone(unittest.TestCase):
	def test1(self):
		print "RUNNING SETUP"
		fn = 'ParticleReceiver1D_standalone.mo'
		sim = simulation.Simulator(fn)
		resfn = sim.model + '_res.mat'
		global RUNSIM
		if RUNSIM=='if-needed':
			RUNSIM = os.path.getctime(resfn) < os.path.getctime(fn)
			print "RUNSIM =",RUNSIM
		if RUNSIM:
			print "COMPILING MODEL"
			sim.compile_model(args=['-d=bltdump'])
			print "COMPILING SIM"
			sim.compile_sim(args=['-s'])
			print "SOLVING MODEL"
			sim.simulate(start=0, stop='1s', step='1s'
				,solver='dassl',nls=None
				#,lv='LOG_DEBUG,LOG_NLS,LOG_SOLVER,LOG_NLS_V'
			)
		self.res = postproc.SimResult(resfn)

		print "RESULTS"

		def getval(n):
			return self.res.interpolate(n,1)

		N = getval('N')
		vl = ['N','T_s[1]','T_s[%d]'%(N,), 'mdot', 'T_in','q_solar','A_ap'
			,'Qdot_rec','Qdot_inc','eta_rec','eps_c[%d]'%(N,)]
		for v in vl:
			print '%s = %f' %(v,getval(v))

		if PLOTME:
			import matplotlib; matplotlib.use('GTKCairo')
			import matplotlib.pyplot as pl
			nr = 1; nc = 6; sp=0
			dx = getval('dx')
			x = [dx*(i+1) for i in range(N)]

			fig = pl.figure()
			sp+=1; ax1 = pl.subplot(nr,nc,sp)
			for vn in ['v_s']:
				Qdotdd = [getval('%s[%d]'%(vn,i+1)) for i in range(N)]
				pl.plot(Qdotdd,x,label=vn,marker='o',markerfacecolor=None)
			pl.legend()
			pl.ylabel('Vertical position / [m]')
			pl.xlabel('Velocity / [m/s]')

			sp+=1; pl.subplot(nr,nc,sp,sharey=ax1)
			for vn in ['th_c']:
				Qdotdd = [getval('%s[%d]'%(vn,i+1))/1000 for i in range(N)]
				pl.plot(Qdotdd,x,label=vn,marker='o',markerfacecolor=None)
			pl.legend()
			pl.xlabel('Thickness / [mm]')

			sp+=1; pl.subplot(nr,nc,sp,sharey=ax1)
			Qdot_new_c = [getval('Qdot_net_c[%d]'%(i+1))/1000. for i in range(N)]
			pl.xlabel('$\dot{Q}_\mathrm{net,c}$ / [kW/m2]')
			pl.plot(Qdot_new_c,x,'ro-',label='backwall')

			sp+=1; pl.subplot(nr,nc,sp,sharey=ax1)
			for vn in ['tau_c']:
				Qdotdd = [getval('%s[%d]'%(vn,i+1)) for i in range(N)]
				pl.plot(Qdotdd,x,label=vn,marker='o',markerfacecolor=None)
			pl.legend()
			pl.xlabel('Fraction / [1]')

			sp+=1; pl.subplot(nr,nc,sp,sharey=ax1)
			for vn in ['eps_c','abs_c']:
				Qdotdd = [getval('%s[%d]'%(vn,i+1)) for i in range(N)]
				pl.plot(Qdotdd,x,label=vn,marker='o',markerfacecolor=None)
			pl.legend()
			pl.xlabel('Fraction / [1]')

			pl.subplot(nr,nc,4,sharey=ax1)
			for vn in ['gc_f','gc_b','g_w','jc_f','jc_b']:
				Qdotdd = [getval('%s[%d]'%(vn,i+1))/1000. for i in range(N)]
				pl.plot(Qdotdd,x,label=vn,marker='o')
			pl.legend()
			pl.xlabel('Heat flux / [kW/m2]')

			sp+=1; pl.subplot(nr,nc,sp,sharey=ax1)
			T_s = [getval('T_s[%d]'%(i+1)) for i in range(N)]
			pl.plot(T_s,x,'ko-',label='particle')
			T_w = [getval('T_w[%d]'%(i+1)) for i in range(N+1)]
			x_w = [getval('x[%d]'%(i+1)) for i in range(N+1)]
			pl.ylim(max(x_w),0)
			pl.plot(T_w,x_w,'ro-',label='backwall')
			pl.xlabel('Temperature / [K]')
			pl.legend()

			pl.show()

		self.assertAlmostEqual(getval('q_solar'),1200*788.8,delta=1.)
		self.assertAlmostEqual(getval('A_ap'),593.7,delta=0.2)
		self.assertAlmostEqual(getval('T_s[1]'),580.3+273.15,delta=2)
		#self.assertAlmostEqual(getval('mdot'),1827,delta=100)
		self.assertAlmostEqual(getval('T_s[%d]'%(N,)),800+273.15,delta=0.5)
		self.assertAlmostEqual(getval('Qdot_rec'),497690.,delta=100.)
		self.assertAlmostEqual(getval('eta_rec'),0.8568,delta=0.0005)

if __name__ == '__main__':
	PLOTME=1
	RUNSIM='if-needed'
	unittest.main()
