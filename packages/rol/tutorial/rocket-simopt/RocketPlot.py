# @HEADER
# *****************************************************************************
#               Rapid Optimization Library (ROL) Package
#
# Copyright 2014 NTESS and the ROL contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

  rocket = open("Rocket.txt")
  lines = [l for l in rocket.readlines()]
  labels = [" ".join(l.split('_')) for l in lines[0].split()]

  data = np.array([[float(x) for x in l.split()] for l in lines[1:]])
  time,bc,hc,bopt,hopt = [data[:,i] for i in range(data.shape[1])]
  dt = time[2]-time[1]

  fig1 = plt.figure(1)
  fig2 = plt.figure(2)

  ax1 = fig1.add_subplot(111)
  ax1.plot(time,bc/dt,label='constant',lw=2) 
  ax1.plot(time,bopt/dt,label='optimal',lw=2) 
  ax1.set_xlabel('time (s)',fontsize=16)
  ax1.set_ylabel('Burn rate (kg/s)', fontsize=16)
  ax1.tick_params(labelsize=14)
  ax1.legend(fontsize=16)
  ax1.set_xlim(time[0],time[-1])
  fig1.savefig("burn_rate.png")

  ax2 = fig2.add_subplot(111)
  ax2.plot(time,hc/1000.0,label='constant',lw=2)
  ax2.plot(time,hopt/1000.0,label='optimal',lw=2)
  ax2.set_xlabel('time (s)',fontsize=16)
  ax2.set_ylabel('altitude (km)', fontsize=16)
  ax2.tick_params(labelsize=14)
  ax2.legend(fontsize=16)
  ax2.set_xlim(time[0],time[-1])
  fig2.savefig("altitude.png")
 
  plt.show()
