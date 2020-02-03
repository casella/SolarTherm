import os
import sys
import time

import numpy as N
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.interpolate import interp1d,interp2d
import matplotlib.cm as cm

from repository.SolsticePy.proces_raw import *
from repository.SolsticePy.cal_layout import radial_stagger
from repository.SolsticePy.cal_field import *
from repository.SolsticePy.cal_sun import *
from repository.SolsticePy.gen_YAML import gen_YAML
from repository.SolsticePy.gen_vtk import *
from repository.SolsticePy.input import Parameters
from repository.SolsticePy.master_crs import *

def set_param(inputs={}):
    '''
    set parameters
    '''

    pm=Parameters()
    for k, v in inputs.iteritems():

        if hasattr(pm, k):
            setattr(pm, k, v)

        else:
            raise RuntimeError("invalid paramter '%s'"%(k,)) 

    pm.dependent_par()
    return pm

def run_simul(inputs={}):
    RAYS=N.r_[1e5]

    pm=set_param(inputs)

    print ''
    print 'Test inputs'
    for k, v in inputs.iteritems():
        print k, '=', getattr(pm, k)
    print ''
    print ''

    TIME=N.array([])

    for r in RAYS:
        print ''

        start=time.time()
        pm.n_rays=int(r)

        casedir=pm.casedir
        pm.saveparam(casedir)

        crs=CRS(casedir)

        crs.annualsolar(nd=int(pm.n_row_oelt), nh=int(pm.n_col_oelt), latitude=float(pm.lat), sunshape=pm.sunshape, sunsize=float(pm.sunsize))

        crs.heliostatfield(field=pm.field_type, num_hst=int(pm.n_helios), hst_w=float(pm.W_helio), hst_h=float(pm.H_helio), hst_z=float(pm.Z_helio), hst_rho=float(pm.rho_helio), slope=float(pm.slope_error), R1=pm.R1, dsep=pm.dsep, tower_h=float(pm.H_tower), tower_r=float(pm.R_tower))

        crs.receiversystem(receiver=pm.rcv_type, rec_w=float(pm.W_rcv), rec_h=float(pm.H_rcv), rec_x=float(pm.X_rcv), rec_y=float(pm.Y_rcv), rec_z=float(pm.Z_rcv), rec_tilt=float(pm.tilt_rcv), rec_grid=int(pm.n_H_rcv), rec_abs=float(pm.alpha_rcv))

        crs.field_design(Q_in_des=pm.Q_in_rcv, latitude=pm.lat, dni_des=pm.dni_des, num_rays=pm.n_rays, genvtk_hst=False)
      

        annualfolder=casedir+'/annual'
        crs.run_annual(annualfolder, num_rays=pm.n_rays, genvtk_hst=False)
                                                                        
        end=time.time()
        print ''
        print ''
        print 'Simulation', int(r) , 'rays'
        print 'total time', end-start, 's' 
        N.savetxt(casedir+'/time.csv', N.r_[r, end-start], fmt='%.4f', delimiter=',')

    tablefile=casedir+'/OELT_Solstice.motab'
    return tablefile

    
    
if __name__=='__main__':
    inputs={'casedir': "./result/demo"}
    run_simul(inputs)


