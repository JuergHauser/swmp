import sys
import os

cwd=os.getcwd()
sys.path.append(cwd+"/../../python")

import swmp

import numpy

wt=swmp.WaveFrontTracker()
wt.read_configuration('input/cur_rat.in')
wt.forward()

wt.read_predictions('output/arrivals.dat')

wt.read_jacobian()
jac=wt.get_jacobian()

vis=swmp.Visualisation()
vis.read_configuration('input/cur_rat.in')
vis.read_raypaths()
vis.read_wavefronts()
fig=vis.get_wavefront_figure(5,5)
fig.savefig('output/wavefronts.png')
