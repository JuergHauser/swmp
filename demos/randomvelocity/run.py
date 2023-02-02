import sys
import os

cwd=os.getcwd()
sys.path.append(cwd+"/../../python")

import swmp

tv=swmp.VelocityModelGenerator()
tv.read_configuration('input/true2dv.in')
tv.run()

sv=swmp.VelocityModelGenerator()
sv.read_configuration('input/start2dv.in')
sv.run()

wt=swmp.WaveFrontTracker()
wt.read_forward_configuration('input/true_rat.in')
wt.forward()
