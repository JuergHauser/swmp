import sys
import os

cwd=os.getcwd()
sys.path.append(cwd+"/../../python")

import swmp

wt=swmp.WaveFrontTracker()
wt.read_configuration('input/current_rat.in')
wt.forward()

vis=swmp.Visualisation()
vis.read_configuration('input/current_rat.in')
vis.read_raypaths()
vis.read_wavefronts()
fig=vis.get_wavefront_figure(5,5)
fig.savefig('output/wavefronts.png')
