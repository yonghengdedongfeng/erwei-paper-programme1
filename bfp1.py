from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core
from pyGDM2 import propagators
from pyGDM2 import linear
from pyGDM2 import tools
from pyGDM2 import visu

import numpy as np
import copy
import matplotlib.pyplot as plt
if __name__ == '__main__':
    mesh = 'hex'
    radius = 100/2. # nm
    step = 10    # nm; a smaller step was used for the paper figure
    geometry = structures.sphere(step, R = radius/step, mesh=mesh)
    geometry.T[0] = geometry.T[0]+100
    material = materials.gold()
    struct = structures.struct(step, geometry, material, normalization=mesh)
    print(" N dipoles:", len(struct.geometry))
    ## --- environment (glass sustrate, vacuum above)
    n1, n2 = 1.0, 1.0
    dyads = propagators.DyadsQuasistatic123(n1=n1, n2=n2)
    NA = 0.90
    f = 1.8     # focal length (mm)
    w0 = f * NA / n2   # n2: structure environment
    wavelengths = [532]
    
    field_generator = fields.HermiteGauss00
    polarization_state= (1,1j,0,0)
    kwargs = dict(xSpot=0.0, ySpot=0.0, zSpot=50.0, kSign = -1.0,
                      NA = NA, f = f, w0 = w0,theta=None, polarization_state= polarization_state, returnField='E')
    
    efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)
    ## --- Simulation object
    sim = core.simulation(struct, efield, dyads)

    sim.scatter()

    Nteta=90
    Nphi=2*72
    teta_below, phi_below, Escat, Etot, E0_below = linear.farfield(
                                    sim, field_index=0,
                                    tetamin=np.pi/2., tetamax=np.pi,
                                    Nteta=Nteta, Nphi=Nphi,return_value="efield")

    tetalist = np.ones((int(Nteta), int(Nphi)))*np.linspace(np.pi/2,np.pi ,int(Nteta))[:,None]
    philist = np.ones((int(Nteta), int(Nphi)))*np.linspace(0, 2*np.pi, int(Nphi), endpoint=False)[None,:]
    Er = (Etot.T[0] * np.sin(tetalist.flatten()) * np.cos(philist.flatten()) + 
                    Etot.T[1] * np.sin(tetalist.flatten()) * np.sin(philist.flatten()) + 
                    Etot.T[2] * np.cos(tetalist.flatten()))
    Ep  = ( Etot.T[0] * np.cos(tetalist.flatten()) * np.cos(philist.flatten()) + 
                    Etot.T[1] * np.sin(philist.flatten()) * np.cos(tetalist.flatten()) - 
                    Etot.T[2] * np.sin(tetalist.flatten()) )
    Es = (-1*Etot.T[0] * np.sin(philist.flatten()) + Etot.T[1] * np.cos(philist.flatten()) )
    
    
    #
    Er = np.abs(Er)**2
    Ep = np.abs(Ep)**2
    Es = np.abs(Es)**2
    #
    #kx/k0 ky/k0
    x1 = np.sin(tetalist.flatten()) * np.cos(philist.flatten())
    y1 = np.sin(tetalist.flatten()) * np.sin(philist.flatten())
   
    delta = 0.025
    num = 81
    x = y = np.linspace(1.0, -1.0, num)
    X, Y = np.meshgrid(x, y)
    Z1 = 0*np.ones((len(x),len(x)),dtype=complex).reshape(X.shape)
    Z2 = 0*np.ones((len(x),len(x)),dtype=complex).reshape(X.shape)
    Z3 = 0*np.ones((len(x),len(x)),dtype=complex).reshape(X.shape)
    for zv in range(len(x)):
        for zh in range(len(y)):
            for index in range(len(x1)):
                #通过对变化范围的减小或增大，我们可以得到不同的分布，这是人为造成的，需要慎重调节
                if np.abs(X[zv,zh]-x1[index])<=0.02 and np.abs(Y[zv,zh]-y1[index])<=0.02:
                    Z1[zv,zh] = E0_below[index,0]
                    Z2[zv,zh] = E0_below[index,1]
                    Z3[zv,zh] = E0_below[index,2]
    Ip = np.abs(Z1)**2 +np.abs(Z2)**2+np.abs(Z3)**2
    fig, ax = plt.subplots(dpi=500)
    im = ax.imshow(Ip,interpolation='bilinear')
    ax.set_title('E0^2')
    ax.set(xlabel='kx/k0', ylabel='ky/k0',title='moment space')
    ax.set_xticks([0,40,81],[-1, 0, 1])
    ax.set_yticks([0,40,80],[1, 0, -1])
    ax.grid()
    fig.colorbar(im, ax=ax, label='')
    plt.show()
