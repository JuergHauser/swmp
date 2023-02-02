
import ctypes

import numpy

class VelocityModel():
    def __init__(self):
        self.x0=None
        self.y0=None
        self.nx=None
        self.ny=None
        self.dx=None
        self.dy=None
        self.cn=None
        self.mv=None
        self.mg=None

    def read(self,fn):
        pass

class TravelTimeData():
    def __init__(self):
        self.obs=None
        self.pred=None
        pass

    def read_observations(self,fn):
        pass

    def read_predictions(self,fn):
        pass

class WaveFrontTracker():

    def __init__(self):
        self.tt=TravelTimeData()
        self.swmp=ctypes.cdll.LoadLibrary("../../modules/libswmp.so")
        pass

    def read_forward_configuration(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.rat_read_conf(fn,ctypes.c_int(len(fn.value)))

    def set_dt(self,dt):
        val = ctypes.c_float(dt)
        self.swmp.set_dt(val)

    def get_dt(self):
        val = ctypes.c_float(-99.9)
        self.swmp.get_dt(ctypes.byref(val))
        return float(val.value)

    def forward(self):
        self.swmp.forward()
        pass

    def read_observations(fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.read_observations(fn,ctypes.c_int(len(fn.value)))
        n = ctypes.c_int(-99)
        self.swmp.get_number_of_observations(ctypes.byref(n))
        tt=numpy.asfortranarray(numpy.zeros([n,6]))
        self.get_observations(numpy_pointer(tt),ctypes.byref(n))
        self.tt.obs=numpy.array(tt)

    def get_data(self):
        return self.tt

    def read_predictions(fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.read_predictions(fn,ctypes.c_int(len(fn.value)))
        n = ctypes.c_int(-99)
        self.swmp.get_number_of_observations(ctypes.byref(n))
        tt=numpy.asfortranarray(numpy.zeros([n,5]))
        self.get_predictions(numpy_pointer(tt),ctypes.byref(n))
        self.tt.pred=numpy.array(tt)

    def get_model_grid(self):
        self.mod=VelocityModel()
        self.swmp.get_model_meta_data(x0,y0,nx,ny,dx,dy,cn)
        self.mod.x0=float(x0.value)
        self.mod.y0=float(y0.value)
        self.mod.nx=int(nx.value)
        self.mod.ny=int(ny.value)
        self.mod.dx=float(dx.value)
        self.mod.dy=float(dy.value)
        self.mod.cn=int(cn.value)
        mg=numpy.asfortranarray(numpy.zeros([(nx+cn*2),(ny+cn*2)]))
        self.get_model_grid(numpy_pointer(mg),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy))
        self.mod.mg=numpy.array(mg)
        return self.mod

    def get_model_vector(self):
        self.mod=VelocityModel()
        self.swmp.get_model_meta_data(x0,y0,nx,ny,dx,dy,cn)
        self.mod.x0=float(x0.value)
        self.mod.y0=float(y0.value)
        self.mod.nx=int(nx.value)
        self.mod.ny=int(ny.value)
        self.mod.dx=float(dx.value)
        self.mod.dy=float(dy.value)
        self.mod.cn=int(cn.value)
        mv=numpy.asfortranarray(numpy.zeros((nx+cn*2)*(ny+cn*2)))
        self.get_model_vector(numpy_pointer(mv),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy))
        self.mod.mv=numpy.array(mv)
        return self.mod

    def set_model_vector(self,mv_):
        mv=numpy.asfortranarray(mv_)
        self.set_model_vector(numpy_pointer(mv),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy))
        return

class VelocityModelGenerator():

    def __init__(self):
        self.swmp=ctypes.cdll.LoadLibrary("../../modules/libswmp.so")
        pass

    def read_configuration(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.gen2dv_read_conf(fn,ctypes.c_int(len(fn.value)))

    def run(self):
        self.swmp.gen2dv_run()

if __name__ == "__main__":
    wt=WavefrontTracker()
    wt.set_dt(1.0)
    val=wt.get_dt()
    print(val)

#fn=ctypes.c_char_p(b'rat.in')
#swmp.read_config_from_file(fn,ctypes.c_int(len(fn.value)))

#dt = ctypes.c_float(0.01)
#swmp.set_dt(ctypes.byref(dt))

#val = ctypes.c_float(0.01)
#swmp.get_dt(val)
#print (val)
