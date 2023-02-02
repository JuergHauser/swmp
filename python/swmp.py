
import ctypes
import numpy

# https://github.com/dhermes/foreign-fortran

def numpy_pointer(array):
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

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
        fn=ctypes.c_char_p((" "*64).encode('UTF-8'))
        self.swmp.get_arrival_prediction_filepath(fn,ctypes.c_int(len(fn.value)))
        fp=(fn.value.decode('UTF-8')).rstrip()
        self.tt.pred=numpy.loadtxt(fp)

        fh=open(fp,'w')
        fh.write("{}\n".format(numpy.shape(self.tt.pred)[0]))
        for i in range(numpy.shape(self.tt.pred)[0]):
            fh.write('{:d} {:d} {:d} {} {} \n'.format(int(self.tt.pred[i,0]),int(self.tt.pred[i,1]),int(self.tt.pred[i,2]),self.tt.pred[i,3],self.tt.pred[i,4]))
        fh.close()

        pass

    def read_observations(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.read_observations(fn,ctypes.c_int(len(fn.value)))
        n = ctypes.c_int(-99)
        self.swmp.get_number_of_observations(ctypes.byref(n))
        tt=numpy.empty([n.value,6], dtype=ctypes.c_float)
        self.swmp.get_observations(tt.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(n))
        self.tt.obs=numpy.array(tt)

    def get_data(self):
        return self.tt

    def read_predictions(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.read_predictions(fn,ctypes.c_int(len(fn.value)))
        n = ctypes.c_int(-99)
        self.swmp.get_number_of_observations(ctypes.byref(n))
        tt=numpy.empty([n.value,5], dtype=ctypes.c_float)
        self.swmp.get_predictions(tt.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(n))
        self.tt.pred=numpy.array(tt)

    def get_model_grid(self):
        self.mod=VelocityModel()
        x0=ctypes.c_float(-99.9)
        y0=ctypes.c_float(-99.9)
        nx=ctypes.c_int(-99)
        ny=ctypes.c_int(-99)
        dx=ctypes.c_float(-99.9)
        dy=ctypes.c_float(-99.9)
        cn=ctypes.c_int(-99)
        self.swmp.get_model_meta_data(x0,y0,nx,ny,dx,dy,cn)
        self.mod.x0=float(x0.value)
        self.mod.y0=float(y0.value)
        self.mod.nx=int(nx.value)
        self.mod.ny=int(ny.value)
        self.mod.dx=float(dx.value)
        self.mod.dy=float(dy.value)
        self.mod.cn=int(cn.value)
        mv=numpy.empty((self.mod.nx+self.mod.cn*2),(self.mod.ny+self.mod.cn*2), dtype=ctypes.c_float)
        self.swmp.get_model_grid(mg.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        self.mod.mg=numpy.array(mg)
        return self.mod.mg

    def get_model_vector(self):
        self.mod=VelocityModel()
        x0=ctypes.c_float(-99.9)
        y0=ctypes.c_float(-99.9)
        nx=ctypes.c_int(-99)
        ny=ctypes.c_int(-99)
        dx=ctypes.c_float(-99.9)
        dy=ctypes.c_float(-99.9)
        cn=ctypes.c_int(-99)
        self.swmp.get_model_meta_data(ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        self.mod.x0=float(x0.value)
        self.mod.y0=float(y0.value)
        self.mod.nx=int(nx.value)
        self.mod.ny=int(ny.value)
        self.mod.dx=float(dx.value)
        self.mod.dy=float(dy.value)
        self.mod.cn=int(cn.value)
        mv=numpy.empty((self.mod.nx+self.mod.cn*2)*(self.mod.ny+self.mod.cn*2), dtype=ctypes.c_float)
        self.swmp.get_model_vector(mv.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        self.mod.mv=numpy.array(mv)
        return self.mod.mv

    def set_model_vector(self,mv_):
        x0=ctypes.c_float(-99.9)
        y0=ctypes.c_float(-99.9)
        nx=ctypes.c_int(-99)
        ny=ctypes.c_int(-99)
        dx=ctypes.c_float(-99.9)
        dy=ctypes.c_float(-99.9)
        cn=ctypes.c_int(-99)
        self.swmp.get_model_meta_data(ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        mv=numpy.asfortranarray(mv_)
        self.swmp.set_model_vector(mv.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
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


class ObservationGenerator():

    def __init__(self):
        self.swmp=ctypes.cdll.LoadLibrary("../../modules/libswmp.so")
        pass

    def read_configuration(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.creobs_read_conf(fn,ctypes.c_int(len(fn.value)))

    def run(self):
        self.swmp.creobs_run()

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
