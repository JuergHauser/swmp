import glob
import os
import os.path

import ctypes

import numpy
import numpy.ctypeslib

import scipy
import scipy.interpolate
import scipy.sparse

import matplotlib
import matplotlib.pyplot
import mpl_toolkits
import mpl_toolkits.basemap


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
        self.m=None


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

    def read_configuration(self,fn_):
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
            fh.write('{} {} {} {} {} \n'.format(int(self.tt.pred[i,0]),int(self.tt.pred[i,1]),int(self.tt.pred[i,2]),self.tt.pred[i,3],self.tt.pred[i,4]))
        fh.close()

        pass

    def read_observations(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.read_observations(fn,ctypes.c_int(len(fn.value)))
        n = ctypes.c_int(-99)
        self.swmp.get_number_of_observations(ctypes.byref(n))
        tt=numpy.empty([n.value,6], dtype=ctypes.c_float)
        self.swmp.get_observations(tt.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(n))
        self.tt.obs=numpy.array(tt,order="C").reshape((6,n.value)).transpose()

    def get_data(self):
        return self.tt

    def read_predictions(self,fn_):
        fn=ctypes.c_char_p(fn_.encode('UTF-8'))
        self.swmp.read_predictions(fn,ctypes.c_int(len(fn.value)))
        n = ctypes.c_int(-99)
        self.swmp.get_number_of_predictions(ctypes.byref(n))
        tt=numpy.empty([n.value,5], dtype=ctypes.c_float)
        self.swmp.get_predictions(tt.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(n))
        self.tt.pred=numpy.array(tt,order="C").reshape((5,n.value)).transpose()

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
        m=numpy.empty((self.mod.nx+self.mod.cn*2)*(self.mod.ny+self.mod.cn*2), dtype=ctypes.c_float)
        self.swmp.get_model_vector(m.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        self.mod.m=numpy.array(m)
        return self.mod.m


    def get_resampled_model_vector(self):
        self.resamod=VelocityModel()
        x0=ctypes.c_float(-99.9)
        y0=ctypes.c_float(-99.9)
        nx=ctypes.c_int(-99)
        ny=ctypes.c_int(-99)
        dx=ctypes.c_float(-99.9)
        dy=ctypes.c_float(-99.9)
        cn=ctypes.c_int(-99)
        self.swmp.get_resampled_model_meta_data(ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        self.resamod.x0=float(x0.value)
        self.resamod.y0=float(y0.value)
        self.resamod.nx=int(nx.value)
        self.resamod.ny=int(ny.value)
        self.resamod.dx=float(dx.value)
        self.resamod.dy=float(dy.value)
        self.resamod.cn=int(cn.value)
        m=numpy.empty((self.resamod.nx+self.resamod.cn*2)*(self.resamod.ny+self.resamod.cn*2), dtype=ctypes.c_float)
        self.swmp.get_resampled_model_vector(m.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        self.resamod.m=numpy.array(m)
        return self.resamod.m

    def set_model_vector(self,m_):
        x0=ctypes.c_float(-99.9)
        y0=ctypes.c_float(-99.9)
        nx=ctypes.c_int(-99)
        ny=ctypes.c_int(-99)
        dx=ctypes.c_float(-99.9)
        dy=ctypes.c_float(-99.9)
        cn=ctypes.c_int(-99)
        self.swmp.get_model_meta_data(ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        m=numpy.asfortranarray(m_)
        self.swmp.set_model_vector(m.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(x0),ctypes.byref(y0),ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(dx),ctypes.byref(dy),ctypes.byref(cn))
        return



    def read_jacobian(self):
        self.swmp.read_jacobian()
        pass

    def get_jacobian(self):
        nr_=ctypes.c_int(-99)
        nc_=ctypes.c_int(-99)
        nnz_=ctypes.c_int(-99)
        self.swmp.get_sparse_jacobian_size(ctypes.byref(nr_),ctypes.byref(nc_),ctypes.byref(nnz_))
        nr=int(nr_.value)
        nc=int(nc_.value)
        nnz=int(nnz_.value)

        jrow_=numpy.empty(nnz, dtype=ctypes.c_float)
        jcol_=numpy.empty(nnz, dtype=ctypes.c_float)
        jval_=numpy.empty(nnz, dtype=ctypes.c_float)

        self.swmp.get_sparse_jacobian(jrow_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),jcol_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),jval_.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),ctypes.byref(nnz_))

        jrow=numpy.array(jrow_)-1
        jcol=numpy.array(jcol_)-1
        jval=numpy.array(jval_)

        return scipy.sparse.csr_array((jval, (jrow, jcol)), shape=(nr, nc))


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


class Visualisation(WaveFrontTracker):

    def __init__(self):
        WaveFrontTracker.__init__(self)
        self.rp=[]
        self.wf=[]
        pass



    def read_raypaths(self):
        fn=ctypes.c_char_p((" "*64).encode('UTF-8'))
        self.swmp.get_raypath_prediction_filepath(fn,ctypes.c_int(len(fn.value)))
        fp=(fn.value.decode('UTF-8')).rstrip()
        fpb=fp.rstrip('.'+os.path.basename(fp).split('.')[-1])
        fps = glob.glob(fpb+'*.'+os.path.basename(fp).split('.')[-1])

        for fp in fps:
            fh=open(fp)
            lines=fh.readlines()
            rp=[]
            irp=[]
            for line in lines:
                fields=line.split()
                if fields[0]==('>') and irp!=[]:
                    irp=numpy.array(irp)
                    rp.append(irp)
                    irp=[]
                elif fields[0]!=('>'):
                    irp.append([float(fields[0]),float(fields[1])])
            if irp!=[]:
                irp=numpy.array(irp)
                rp.append(irp)

            self.rp.append(rp)

    def read_wavefronts(self):
        fn=ctypes.c_char_p((" "*64).encode('UTF-8'))
        self.swmp.get_wavefront_prediction_filepath(fn,ctypes.c_int(len(fn.value)))
        fp=(fn.value.decode('UTF-8')).rstrip()
        fh=open(fp)
        lines=fh.readlines()
        wf=[]
        for line in lines:
            fields=line.split()
            if fields[0]==('>') and wf!=[]:
                wf=numpy.array(wf)
                self.wf.append(wf)
                wf=[]
            elif fields[0]!=('>'):
                wf.append([float(fields[0]),float(fields[1])])
        if  wf!=[]:
            wf=numpy.array(wf)
            self.wf.append(wf)

    def get_model_figure(self,nrx_,nry_):
        nrx=ctypes.c_int(nrx_)
        nry=ctypes.c_int(nry_)
        self.swmp.resample_model(nrx,nry)
        m = self.get_resampled_model_vector()
        x0=self.resamod.x0
        y0=self.resamod.y0
        nx=self.resamod.nx
        ny=self.resamod.ny
        x1=x0+self.resamod.dx*self.resamod.nx
        y1=y0+self.resamod.dy*self.resamod.ny
        xc=(x0+x1)/2.0
        yc=(y0+y1)/2.0

        x = numpy.linspace(x0, x1, nx)
        y = numpy.linspace(y0, y1, ny)

        yy, xx = numpy.meshgrid(y, x)
        zz =numpy.reshape(m,(nx,ny))

        graph = mpl_toolkits.basemap.Basemap (llcrnrlon=x0,llcrnrlat=y0,urcrnrlon=x1,urcrnrlat=y1,
            resolution='i',projection='merc',lon_0=xc,lat_0=yc)

        cmap = matplotlib.pyplot.colormaps['Greys_r']
        graph.pcolormesh(xx, yy, zz,latlon=True,cmap=cmap)
        graph.drawcoastlines()

        parallels = numpy.arange(-81,81,10)
        graph.drawparallels(parallels,labels=[True,False,False,True])
        meridians = numpy.arange(10,351,20)
        graph.drawmeridians(meridians,labels=[False,True,True,False])

        return matplotlib.pyplot.gcf()


    def get_raypath_figure(self,nrx_,nry_):
        nrx=ctypes.c_int(nrx_)
        nry=ctypes.c_int(nry_)
        self.swmp.resample_model(nrx,nry)
        m = self.get_resampled_model_vector()
        x0=self.resamod.x0
        y0=self.resamod.y0
        nx=self.resamod.nx
        ny=self.resamod.ny
        x1=x0+self.resamod.dx*self.resamod.nx
        y1=y0+self.resamod.dy*self.resamod.ny
        xc=(x0+x1)/2.0
        yc=(y0+y1)/2.0

        x = numpy.linspace(x0, x1, nx)
        y = numpy.linspace(y0, y1, ny)

        yy, xx = numpy.meshgrid(y, x)
        zz =numpy.reshape(m,(nx,ny))

        graph = mpl_toolkits.basemap.Basemap (llcrnrlon=x0,llcrnrlat=y0,urcrnrlon=x1,urcrnrlat=y1,
            resolution='i',projection='merc',lon_0=xc,lat_0=yc)

        cmap = matplotlib.pyplot.colormaps['Greys_r']
        graph.pcolormesh(xx, yy, zz,latlon=True,cmap=cmap)
        graph.drawcoastlines()
        # http://tsitsul.in/blog/coloropt/
        colors=['#ebac23','#b80058','#008cf9','#006e00','#00bbad','#d163e6','#b24502','#ff9287','#5954d6','#00c6f8','#878500','#00a76c']
        ic=0
        for rp in self.rp:
            for irp in rp:
                graph.plot(irp[:,0],irp[:,1],linewidth=1,color=colors[ic],latlon=True)
            ic=ic+1

        parallels = numpy.arange(-81,81,10)
        graph.drawparallels(parallels,labels=[True,False,False,True])
        meridians = numpy.arange(10,351,20)
        graph.drawmeridians(meridians,labels=[False,True,True,False])

        return matplotlib.pyplot.gcf()

    def get_wavefront_figure(self,nrx_,nry_):
        nrx=ctypes.c_int(nrx_)
        nry=ctypes.c_int(nry_)
        self.swmp.resample_model(nrx,nry)
        m = self.get_resampled_model_vector()
        x0=self.resamod.x0
        y0=self.resamod.y0
        nx=self.resamod.nx
        ny=self.resamod.ny
        x1=x0+self.resamod.dx*self.resamod.nx
        y1=y0+self.resamod.dy*self.resamod.ny
        xc=(x0+x1)/2.0
        yc=(y0+y1)/2.0

        x = numpy.linspace(x0, x1, nx)
        y = numpy.linspace(y0, y1, ny)

        yy, xx = numpy.meshgrid(y, x)
        zz =numpy.reshape(m,(nx,ny))

        graph = mpl_toolkits.basemap.Basemap (llcrnrlon=x0,llcrnrlat=y0,urcrnrlon=x1,urcrnrlat=y1,
            resolution='i',projection='merc',lon_0=xc,lat_0=yc)

        cmap = matplotlib.pyplot.colormaps['Greys_r']
        graph.pcolormesh(xx, yy, zz,latlon=True,cmap=cmap)
        graph.drawcoastlines()
        # http://tsitsul.in/blog/coloropt/
        colors=['#ebac23','#b80058','#008cf9','#006e00','#00bbad','#d163e6','#b24502','#ff9287','#5954d6','#00c6f8','#878500','#00a76c']
        ic=0

        for wf in self.wf:
            graph.plot(wf[:,0],wf[:,1],linewidth=1,color='k',latlon=True)

        for rp in self.rp:
            for irp in rp:
                graph.plot(irp[:,0],irp[:,1],linewidth=1,color=colors[ic],latlon=True)
            ic=ic+1

        parallels = numpy.arange(-81,81,10)
        graph.drawparallels(parallels,labels=[True,False,False,True])
        meridians = numpy.arange(10,351,20)
        graph.drawmeridians(meridians,labels=[False,True,True,False])

        return graph



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
