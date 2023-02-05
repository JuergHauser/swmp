

import os
import sys
import subprocess

# We check if swmp is actually available otherwise we pull and build it.
def build_swmp():
    path = os.path.dirname(__file__)
    proc = subprocess.Popen("cmake .\n", cwd=path, shell=True)
    proc.wait()
    proc = subprocess.Popen("make \n", cwd=path, shell=True)
    proc.wait()

path = os.path.dirname(__file__)
sys.path.append(path+"/_deps/swmp-build/python")

if not(os.path.exists(path+"/_deps/swmp-build/python/swmp.py")):
    build_swmp()

import swmp

from cofi_espresso import EspressoProblem
from cofi_espresso.exceptions import InvalidExampleError

class SurfaceWaveMultipathing(EspressoProblem):
    """Forward simulation class
    """

    # TODO fill in the following metadata.
    metadata = {
        "problem_title": "SurfaceWaveMultipathing",                # To be u
        "problem_short_description": "",    # 1-3 sentences

        "author_names": ["Juerg Hauser"],    # List of names e.g. author_names = ["Sally Smith", "Mark Brown"]

        "contact_name": "Juerg Hauser",         # Contact for contributor/maintainer of espresso example
        "contact_email": "juerg.hauser@csiro.au",

        "citations": [("Hauser, J., Sambridge, M. and Rawlinson, N. (2008). Multiarrival wavefront tracking and its applications. Geochem. Geophys. Geosyst., 9(11), Q11001. https://doi.org/10.1029/2008GC002069","")], # Reference to publication(s) that describe this example. In most
                                # cases there will only be a single entry in this list.
                                # List of (citation, doi) pairs e.g.
                                # citations = [("Newton, I (1687). Philosophiae naturalis principia mathematica.", "")]
                                # If there are no citations, use empty list `[]`

        "linked_sites": [("Parent project on github","https://github.com/JuergHauser/swmp")],  # List of (title, address) pairs for any websites that
                                    # should be linked in the documentation, e.g.
                                    # linked_sites = [("Parent project on Github","https://github.com/user/repo"),
                                    #                 ("Data source"),"https://www.data.com") ]
                                    # If there are no links, use empty list `[]`
    }


    def __init__(self, example_number=1):
        super().__init__(example_number)

        """you might want to set some useful example-specific parameters here
        """
        if example_number == 1:
        # swmp random velocity model example
            self.swmp_demo='randomvelocity'
            path = os.path.dirname(__file__)
            self.wdir = path+'/_deps/swmp-build/demos/randomvelocity'
            os.chdir(self.wdir)

            tv=swmp.VelocityModelGenerator()
            tv.read_configuration('input/true2dv.in')
            tv.run()

            sv=swmp.VelocityModelGenerator()
            sv.read_configuration('input/start2dv.in')
            sv.run()

            cv=swmp.VelocityModelGenerator()
            cv.read_configuration('input/current2dv.in')
            cv.run()


            self.true=swmp.WaveFrontTracker()
            self.true.read_forward_configuration('input/current_rat.in')

            self.start=swmp.WaveFrontTracker()
            self.start.read_forward_configuration('input/current_rat.in')

            self.current=swmp.WaveFrontTracker()
            self.current.read_forward_configuration('input/current_rat.in')


            #  check if the data file exists if not generate it
            if not(os.path.exists(path+'/_deps/swmp-build/demos/randomvelocity/output/data/observed.dat')):
                self.true.forward()
                og=swmp.ObservationGenerator()
                og.read_configuration('input/creobs.in')
                og.run()

        # else:
        #     raise InvalidExampleError



    @property
    def description(self):
        raise NotImplementedError               # optional

    @property
    def model_size(self):
        return len(self.current.get_model_vector())

    @property
    def data_size(self):
        return len(self.data())

    @property
    def good_model(self):
        return self.true.get_model.vector()

    @property
    def starting_model(self):
        return self.start.get_model.vector()

    @property
    def data(self):
        raise NotImplementedError               # TODO implement me
        self.true.read_observations(self.wdir+'/output/data/observations')
        return self.true.tt.obs[:,4]


    @property
    def covariance_matrix(self):                # optional
        raise NotImplementedError

    @property
    def inverse_covariance_matrix(self):
        raise NotImplementedError               # optional

    def forward(self, model, with_jacobian=False):
        if with_jacobian:
            raise NotImplementedError           # optional
        else:
            self.current.forward()
            self.current.read_predictions(self.wdir+'/output/current/predictions')
            return self.tt.current.preds[:,4]

    def jacobian(self, model):
        raise NotImplementedError               # optional

    def plot_model(self, model):
        raise NotImplementedError               # optional

    def plot_data(self, data, data2=None):
        raise NotImplementedError               # optional

    def misfit(self, data, data2):
        raise NotImplementedError               # optional

    def log_likelihood(self,data,data2):
        raise NotImplementedError               # optional

    def log_prior(self, model):
        raise NotImplementedError               # optional

