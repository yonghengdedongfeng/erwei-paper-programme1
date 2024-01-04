"""
Created on Thu May 25 13:47:58 2023
这个文件是理论计算焦场的最新进度文件，计算的是圆偏光焦场
@author: Administrator
"""

from pyGDM2 import structures
from pyGDM2 import visu
import numpy as np
from pyGDM2 import tools
from pyGDM2 import propagators
from pyGDM2 import fields
from pyGDM2 import linear
from pyGDM2 import core
from pyGDM2 import materials
import matplotlib.pyplot as plt
import warnings
import types
import copy
#import radial_pol_beams
if __name__ == '__main__':
    #%%
    def _automatic_projection_select(X,Y,Z):
        if  len(np.unique(Z))==1:
            projection = 'xy'
        elif  len(np.unique(Y))==1:
            projection = 'xz'
        elif  len(np.unique(X))==1:
            projection = 'yz'
        else:
            global WARNING_AUTO_PROJECTION_OFF 
            if not WARNING_AUTO_PROJECTION_OFF:
                warnings.warn("3D data. Falling back to XY projection...")
            projection = 'xy'   # fallback
            WARNING_AUTO_PROJECTION_OFF = True
        
        return projection
## 2D colorplot
    def vectorfield_color(NF, struct=None, projection='auto', complex_part='real', 
                          tit='', ax=None, 
                          slice_level=-9999, NX=None, NY=None, fieldComp='I', 
                          borders=0, clim=None, grid_interpolation='linear',
                          show=True, return_map_array=False,
                          **kwargs):
    
        if ax is None:
            ax = plt.gca()
        
        NF_plt = copy.deepcopy(NF)
        if len(NF) == 2:
            NF_plt = NF_plt[1]
        
        if len(NF_plt.T) == 6:
            X,Y,Z, Ex,Ey,Ez = np.transpose(NF_plt)
        elif len(NF_plt.T) == 3 and struct is not None:
            Ex,Ey,Ez = np.transpose(NF_plt)
            X,Y,Z = tools.get_geometry(struct)
        else:
            raise ValueError("Error: Wrong number of columns in vector field. Expected (Ex,Ey,Ez)-tuples + `simulation` object or (x,y,z, Ex,Ey,Ez)-tuples.")
        
        X = np.real(X)
        Y = np.real(Y)
        Z = np.real(Z)
        step = tools.get_step_from_geometry(np.transpose([X,Y,Z]))
        
        if projection.lower() == 'auto':
            projection = _automatic_projection_select(X,Y,Z)
        
        if fieldComp.lower() not in ['i', 'abs']:
            if complex_part.lower() == "real":
                Ex, Ey, Ez = Ex.real, Ey.real, Ez.real
            elif complex_part.lower() == "imag":
                Ex, Ey, Ez = Ex.imag, Ey.imag, Ez.imag
            elif complex_part.lower() == "compabs":
                Ex, Ey, Ez = np.abs(Ex)/np.abs(Ey), np.abs(Ey), np.abs(Ez)
            elif complex_part.lower() == "phi":
                Ex, Ey, Ez = 2*np.imag(Ex*Ey.conjugate())/(np.abs(Ey)**2+np.abs(Ex)**2), np.angle(Ey), np.angle(Ez)
                
            else:
                raise ValueError("Error: Unknown `complex_part` argument value." +
                                 " Must be either 'real' or 'imag'.")
        
        ######后续修改过添加了‘iex','iey','iez'
        if fieldComp.lower() in ['i', 'abs']:
            EF = np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2
            if fieldComp.lower() == 'abs':
                EF = np.sqrt(EF)
        elif fieldComp.lower() == 'ex':
            EF = Ex
        elif fieldComp.lower() == 'ey':
            EF = Ey
        elif fieldComp.lower() == 'ez':
            EF = Ez
        elif fieldComp.lower() == 'cex':
            EF = np.abs(Ex)**2
        elif fieldComp.lower() == 'cey':
            EF = np.abs(Ey)**2
        elif fieldComp.lower() == 'cez':
            EF = np.abs(Ez)**2
        elif fieldComp.lower() == 'phix':
            EF = Ex
        elif fieldComp.lower() == 'phiy':
            EF = Ey
        elif fieldComp.lower() == 'phiz':
            EF = Ez
          
        
        ## map / slicing
        if projection.lower() == 'xy':
            SLICE = copy.deepcopy(Z)
        if projection.lower() == 'xz':
            SLICE = copy.deepcopy(Y)
        if projection.lower() == 'yz':
            SLICE = copy.deepcopy(X)
            
        slice_level_closest = np.unique(SLICE)[np.argmin(np.abs(np.unique(SLICE) - slice_level))]
        if slice_level_closest != slice_level:
            if np.abs(slice_level) != 9999:
                warnings.warn("slice level Z={} contains no meshpoints. Using closest level containing meshpoint at Z={}".format(
                                              slice_level, slice_level_closest))
            slice_level = slice_level_closest
        if projection.lower() == 'xy':
            XYZList = np.transpose([X[(SLICE == slice_level)],Y[(SLICE == slice_level)], EF[(SLICE == slice_level)]])
        if projection.lower() == 'xz':
            XYZList = np.transpose([X[(SLICE == slice_level)],Z[(SLICE == slice_level)], EF[(SLICE == slice_level)]])
        if projection.lower() == 'yz':
            XYZList = np.transpose([Y[(SLICE == slice_level)],Z[(SLICE == slice_level)], EF[(SLICE == slice_level)]])
        
        
        MAP, extent = tools.list_to_grid(XYZList, NX, NY, interpolation=grid_interpolation)
        
        if tit: ax.set_title(tit)
        img = ax.imshow(MAP, extent=extent, **kwargs)
        if clim: img.set_clim(clim)
        
        ax.set_xlim(extent[0] - borders*step, extent[1] + borders*step)
        ax.set_ylim(extent[2] - borders*step, extent[3] + borders*step)
        
        if show: 
            if fieldComp.lower() == 'i':
                plt.colorbar(img, label=r'$|E|^2 / |E_0|^2$')
            else:
                plt.colorbar(img, label=r'$E / |E_0|$')
            ax.set_aspect("equal")
            ax.set_xlabel("{} (nm)".format(projection[0]))
            ax.set_ylabel("{} (nm)".format(projection[1]))
            plt.show()
        else:
            if return_map_array:
                return img, MAP
            else:
                return img
    
    Npt = 21 
    #产生坐标
    r_probe_XY = tools.generate_NF_map(-100,100,Npt, -100,100,Npt,  Z0=80, projection='XY')
    wavelength = 633
    n1 = 1.0
    n2 = n3 = 1.0    # homogeneous
    spacing = 100000.0      # thickness of layer n2 (in nm)
    env_dict = dict(eps1=n1**2, eps2=n2**2, eps3=n3**2, spacing=spacing)
    
    NA = 0.9  # focusing numerical aperture
    f = 1.8   # lens focal distance (mm)
    f0 = 1     # filling factor
    w0 = f0*f*NA/n2
    #####圆偏光定义
    ####left hand rotating light
    polarization_state= (1,1j,0,0)
    kwargs = dict(xSpot=0.0, ySpot=0.0, zSpot=80.0, kSign = -1.0,
                  NA = NA, f = f, w0 = w0,theta=None, polarization_state= polarization_state, returnField='E')
    field_generator = fields.focused_beams.HermiteGauss00
    NF_XY = field_generator(r_probe_XY, env_dict, wavelength,**kwargs)
    x = copy.deepcopy(NF_XY)
  
    #第二次归一化
    Ex,Ey,Ez = np.transpose(x)
    mochang = np.emath.sqrt(np.conjugate(Ex)*Ex+np.conjugate(Ey)*Ey+np.conjugate(Ez)*Ez)
    Edot  = np.emath.sqrt(np.square(Ex)+np.square(Ey)+np.square(Ez))
    phig = -1*np.angle(Edot)/2
    Ex = Ex.reshape(shap2.shape)*np.exp(-1j*phig).reshape(shap2.shape)/mochang.reshape(shap2.shape)
    Ey = Ey.reshape(shap2.shape)*np.exp(-1j*phig).reshape(shap2.shape)/mochang.reshape(shap2.shape)
    Ez = Ez.reshape(shap2.shape)*np.exp(-1j*phig).reshape(shap2.shape)/mochang.reshape(shap2.shape)
    
    #画偏振椭圆图
    fig = plt.figure(dpi = 500)
    ax = fig.add_subplot(111, projection='3d')
    # 使用plot_surface方法绘制椭圆
    for i in range(441):
        wt = np.linspace(0, 2 * np.pi, 100)
        def ellipse(wt):
            #这里的系数只是为了放大一下椭圆大小，观感好看
            x = 5*np.real(Ex[0 + np.divmod(i,21)[0],0 + np.divmod(i,21)[1]] * np.exp(-1j*wt))
            y = 5*np.real(Ey[0 + np.divmod(i,21)[0],0 + np.divmod(i,21)[1]] * np.exp(-1j*wt))
            z = 5*np.real(Ez[0 + np.divmod(i,21)[0],0 + np.divmod(i,21)[1]] * np.exp(-1j*wt))
            return x, y, z
        x, y, z = ellipse(wt)
        ax.plot(x-100+10*np.divmod(i,21)[0], y-100+10*np.divmod(i,21)[1], z, color='r',linewidth=0.5)
    
    # 设置坐标轴的范围，使图形完整显示
    ax.set_xlim([-120,120])
    ax.set_ylim([-120, 120])
    ax.set_zlim([-10, 10])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(True)
    # 显示图形
    plt.show()



  #其他的画图方式
  projection='XY'
  Nskipvec = 1
  plt.figure(figsize=(5, 5),dpi = 500)
  plt.subplot(111, aspect='equal')
  v = vectorfield_color(NF_XY, r_probe_XY, complex_part="imag",projection=projection, fieldComp='ex', tit='XY-Plane-real(Ex)', show=0)
  plt.colorbar(v)
  plt.figure(figsize=(5, 5))
  plt.subplot(111, aspect='equal')
  v = vectorfield_color(NF_XY, r_probe_XY, complex_part="imag",projection=projection, fieldComp='cey', tit='XY-Plane-imag(Ey)', show=0)
  plt.colorbar(v)
  plt.figure(figsize=(5, 5),dpi = 500)
  plt.subplot(111, aspect='equal')
  v = vectorfield_color(NF_XY, r_probe_XY, complex_part="imag",projection=projection, fieldComp='ez', tit='XY-Plane-abs(Ez)', show=0)
  plt.colorbar(v)
  plt.figure(figsize=(5, 5))
  plt.subplot(111, aspect='equal')
  v = vectorfield_color(NF_XY, r_probe_XY, projection=projection, fieldComp='i', tit='XY-Plane-I', show=0)
  plt.colorbar(v)
