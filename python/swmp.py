
import ctypes
import pathlib

class WaveFrontTracker():

    def __init__(self,cfn=None):
        self.swmp=ctypes.cdll.LoadLibrary("../../modules/libswmp.so")
        if cfn!=None:
            fn=ctypes.c_char_p(cfn.encode('UTF-8'))
            self.swmp.rat_read_conf(fn,ctypes.c_int(len(fn.value)))

    def set_dt(self,dt):
        val = ctypes.c_float(dt)
        self.swmp.set_dt(val)

    def get_dt(self):
        val = ctypes.c_float(-99.9)
        self.swmp.get_dt(ctypes.byref(val))
        return float(val.value)

    def run(self):
        self.swmp.rat_run()
        pass


class VelocityModelGenerator():

    def __init__(self,cfn=None):
        self.swmp=ctypes.cdll.LoadLibrary("../../modules/libswmp.so")
        if cfn!=None:
            fn=ctypes.c_char_p(cfn.encode('UTF-8'))
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
