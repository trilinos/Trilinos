# @HEADER
# ************************************************************************
#
#               Rapid Optimization Library (ROL) Package
#                 Copyright (2014) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact lead developers:
#              Drew Kouri   (dpkouri@sandia.gov) and
#              Denis Ridzal (dridzal@sandia.gov)
#
# ************************************************************************
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
