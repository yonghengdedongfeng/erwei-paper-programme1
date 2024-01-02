# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 16:46:11 2023
这段程序把球坐标画图改为了笛卡尔坐标系下的图像
主要用kx,ky动量空间表示后焦面
@author: Administrator
"""


import numpy as np

import matplotlib.pyplot as plt


if __name__ == '__main__':
    Nteta=90
    Nphi=2*72
    tetalist = np.ones((int(Nteta), int(Nphi)))*np.linspace(np.pi/2,np.pi ,int(Nteta))[:,None]
    philist = np.ones((int(Nteta), int(Nphi)))*np.linspace(0, 2*np.pi, int(Nphi), endpoint=False)[None,:]
    E0_below = np.arange(90 * 144)
    
    
    x1 = np.sin(tetalist.flatten()) * np.cos(philist.flatten())
    y1 = np.sin(tetalist.flatten()) * np.sin(philist.flatten())
    delta = 0.025
    num = 81
    x = y = np.linspace(1.0, -1.0, 81)
    X, Y = np.meshgrid(x, y)
    Z = 0*np.ones((len(x),len(x)),dtype=complex).reshape(X.shape)
    for zv in range(len(x)):
        for zh in range(len(y)):
            for index in range(len(x1)):
                if np.abs(X[zv,zh]-x1[index])<=delta and np.abs(Y[zv,zh]-y1[index])<=delta:
                    Z[zv,zh] = E0_below[index]
    Ip = np.abs(Z)**2
    fig, ax = plt.subplots()
    im = ax.imshow(Ip)
    plt.show()