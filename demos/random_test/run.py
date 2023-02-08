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
wt.read_configuration('input/true_rat.in')
wt.forward()

og=swmp.ObservationGenerator()
og.read_configuration('input/creobs.in')
og.run()


wt.read_predictions('output/true/arrivals.dat')
wt.read_observations('output/data/observed.dat')

tt=wt.get_data()

print("Observations")
print(tt.obs)
print("Predictions")
print(tt.pred)
print("Jacobian")
wt.read_jacobian()
jac=wt.get_jacobian()
print(jac)
#vis=swmp.Visualisation()
#vis.read_configuration('input/true_rat.in')
#vis.read_raypaths()
#fig=vis.get_raypath_figure(5,5)
#fig.savefig('output/true/raypaths.png')
