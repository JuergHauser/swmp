import sys
import os

cwd=os.getcwd()
sys.path.append(cwd+"/../../python")

import swmp

tv=swmp.VelocityModelGenerator()
tv.read_configuration('input/true_2dv.in')
tv.run()

sv=swmp.VelocityModelGenerator()
sv.read_configuration('input/start_2dv.in')
sv.run()

sv=swmp.VelocityModelGenerator()
sv.read_configuration('input/current_2dv.in')
sv.run()

wt=swmp.WaveFrontTracker()
wt.read_configuration('input/true_rat.in')
wt.forward()

vis=swmp.Visualisation()
vis.read_configuration('input/true_rat.in')
vis.read_raypaths()
fig=vis.get_raypath_figure(5,5)
fig.savefig('output/true/raypaths.png')

og=swmp.ObservationGenerator()
og.read_configuration('input/creobs.in')
og.run()

wt=swmp.WaveFrontTracker()
wt.read_configuration('input/current_rat.in')
wt.read_observations('data/observed.dat')
wt.forward()
wt.read_predictions('output/current/arrivals.dat')
wt.read_jacobian()
wt.join_observations_and_predictions()

print("observations")
print(wt.get_observations())

print("predictions")
print(wt.get_predictions())

print("joint observations")
print(wt.get_joint_observations())

print("joint predictions")
print(wt.get_joint_predictions())

print("jacobian")
print(wt.get_jacobian())

print("joint jacobian")
print(wt.get_joint_jacobian())

wt=swmp.WaveFrontTracker()


