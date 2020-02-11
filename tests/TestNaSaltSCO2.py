#! /bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
import unittest

from solartherm import simulation
from solartherm import postproc

from math import pi

import os

class TestScheduler(unittest.TestCase):
	def setUp(self):
		fn = '../examples/NaSaltsCO2System.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model()
		sim.compile_sim(args=['-s'])
		sim.simulate(start=0, stop='1y', step='120s',solver='dassl')
		self.res = postproc.SimResultElec(sim.res_fn)
		self.perf = self.res.calc_perf()

	def test_sched(self):
		def getval(n):
			return self.res.interpolate(n,1)

		print "-----------------------------------------------------------------------"
		print 'Energy per year:                          %6.2f MWh'%(self.perf[0])
		print 'Capacity factor:                          %6.2f %%'%(self.perf[2])
		print 'LCOE:                                     $%6.2f/MWh'%(self.perf[1])
		print "-----------------------------------------------------------------------"
		print 'Receiver thermal input at design point:   %6.2f MW'%(getval('R_des')/1e6)
		print 'Receiver thermal output at design point:  %6.2f MW'%(getval('Q_rec_out')/1e6)
		print 'Power block gross rating at design point: %6.2f MW'%(getval('P_gross')/1e6)
		print 'Solar multiple:                           %4.1f   '%(getval('SM'))
		print 'Receiver diameter:                        %4.2f  m'%(getval('D_receiver'))
		print 'Receiver height:                          %4.2f  m'%(getval('H_receiver'))
		print 'Tower height:                             %4.1f  m'%(getval('H_tower'))
		print 'Number of heliostats:                     %4.1f  m'%(getval('n_heliostat'))

		f = open('report.tex','w+')
		content = r'''\documentclass{article}
\begin{document}

\section{Sodium-Salt single tower system}

\subsection{High-level performance parameters}

\begin{center}
\begin{tabular}{ ll } 
	\hline\hline
	Energy per year (MWh):	& %6.2f \\
	Capacity factor (\%%):	& %6.2f \\
	LCOE (\$/MWh):		& %6.2f \\
	\hline\hline
\end{tabular}
\end{center}

\subsection{High-level design parameters}

\begin{center}
\begin{tabular}{ ll } 
	\hline\hline
	Receiver thermal input at design point (MWt):		& %6.1f \\
	Receiver thermal output at design point (MWt):		& %6.1f \\
	Power block gross rating at design point (MWe):		& %6.1f \\
	Full load hours of storage (h):				& %6.1f \\
	Solar multiple:						& %6.1f \\
	Receiver diameter (m):					& %6.1f \\
	Receiver height (m):					& %6.1f \\
	Tower height (m):					& %6.1f \\
	Number of heliostats:					& %i \\
	\hline\hline
\end{tabular}
\end{center}

\subsection{Costing assuptions}

\begin{center}
\begin{tabular}{ ll } 
	\hline\hline
	Real discount rate (\%%):					& %2.1f \\
	Inflation rate rate (\%%):					& %2.1f \\
	Lifetime of plant (years):					& %i \\
	Construction time (years):					& %i \\
	Subsidies on initial investment costs (\$):			& %6.1f \\
	Field cost per design aperture area (\$/m$^2$):			& %6.1f \\
	Site improvements cost per area (\$/m$^2$):			& %6.1f \\
	Storage cost per energy capacity (\$/kWht):			& %6.1f \\
	Power block cost per gross rated power (\$/kWe):		& %6.1f \\
	Balance of plant cost per gross rated power (\$/kWe):		& %6.1f \\
	Land cost per area (\$/acre):					& %6.1f \\
	Fixed O\&M cost per nameplate per year (\$/year):		& %6.1f \\
	Variable O\&M cost per production per year (\$/kWe/year):	& %6.1f \\
	\hline\hline
\end{tabular}
\end{center}

\subsection{Equipment and maintenance costs}

\begin{center}
\begin{tabular}{ ll } 
	\hline\hline
	Field cost (M\$):			& %6.1f \\
	Site improvements cost (M\$):		& %6.1f \\
	Tower cost (M\$):			& %6.1f \\
	Receiver cost (M\$):			& %6.1f \\
	Heat exchanger cost (M\$):		& %6.1f \\
	Storage cost (M\$):			& %6.1f \\
	Power block cost (M\$):			& %6.1f \\
	Balance of plant cost (M\$):		& %6.1f \\
	Piping cost (M\$):			& %6.1f \\
	Cold Pumps cost (M\$):			& %6.1f \\
	Direct capital cost subtotal (M\$):	& %6.1f \\
	Contingency cost (M\$):			& %6.1f \\
	Total direct capital cost (M\$):	& %6.1f \\
	EPC and owner costs (M\$):		& %6.1f \\
	Land cost (M\$):			& %6.1f \\
	Total capital (installed cost) (M\$):	& %6.1f \\
	\hline\hline
\end{tabular}
\end{center}

\end{document}'''

		f.write(content%(
		self.perf[0],
		self.perf[2],
		self.perf[1],
		getval('R_des')/1e6,
		getval('Q_rec_out')/1e6,
		getval('P_gross')/1e6,
		getval('t_storage'),
		getval('SM'),
		getval('D_receiver'),
		getval('H_receiver'),
		getval('H_tower'),
		getval('n_heliostat'),
		getval('r_disc')*100,
		getval('r_i')*100,
		getval('t_life'),
		getval('t_cons'),
		getval('f_Subs'),
		getval('pri_field'),
		getval('pri_site'),
		getval('pri_storage')*1e3*3600,
		getval('pri_block')*1e3,
		getval('pri_bop')*1e3,
		getval('pri_land')*4046.86,
		getval('pri_om_name')*1e3,
		getval('pri_om_prod')*1e6*3600,
		getval('C_field')/1e6,
		getval('C_site')/1e6,
		getval('C_tower')/1e6,
		getval('C_receiver')/1e6,
		getval('C_hx')/1e6,
		getval('C_storage')/1e6,
		getval('C_block')/1e6,
		getval('C_bop')/1e6,
		getval('C_piping')/1e6,
		getval('C_pumps')/1e6,
		getval('C_cap_dir_sub')/1e6,
		getval('C_contingency')/1e6,
		getval('C_cap_dir_tot')/1e6,
		getval('C_EPC')/1e6,
		getval('C_land')/1e6,
		getval('C_cap')/1e6))
		f.close()

if __name__ == '__main__':
	unittest.main()
