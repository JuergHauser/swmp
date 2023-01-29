import sys
import os

cwd=os.getcwd()
sys.path.append(cwd+"/../../python")

import swmp

tv=swmp.VelocityModelGenerator(cfn='input/true2dv.in')
tv.run()

sv=swmp.VelocityModelGenerator(cfn='input/start2dv.in')
sv.run()

wt=swmp.WaveFrontTracker(cfn='input/true_rat.in')
wt.run()
