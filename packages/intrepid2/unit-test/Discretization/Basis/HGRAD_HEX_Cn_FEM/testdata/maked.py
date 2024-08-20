# @HEADER
# *****************************************************************************
#                           Intrepid2 Package
#
# Copyright 2007 NTESS and the Intrepid2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER
psi0 = lambda x: 0.5 * x * ( x - 1.0 )
psi1 = lambda x: ( 1.0 - x ) * ( 1.0 + x )
psi2 = lambda x: 0.5 * x * ( 1.0 + x )
psi = [psi0,psi1,psi2]

dpsi0 = lambda x: x - 0.5
dpsi1 = lambda x: -2.0 * x
dpsi2 = lambda x: x + 0.5
dpsi = [dpsi0,dpsi1,dpsi2]

ddpsi0 = lambda x: 1.0
ddpsi1 = lambda x: -2.0
ddpsi2 = lambda x: 1.0
ddpsi = [ddpsi0, ddpsi1, ddpsi2]

xs = [ -1.0 , 0.0 , 1.0 ]
pts = [ (xs[i],xs[j],xs[k]) for k in range(3) for j in range(3) for i in range(3) ]


# This is D1
#for k in range(3):
#    for j in range(3):
#        for i in range(3):
#            for pt in pts:
#                print dpsi[i](pt[0])*psi[j](pt[1])*psi[k](pt[2])
#                print psi[i](pt[0])*dpsi[j](pt[1])*psi[k](pt[2])
#                print psi[i](pt[0])*psi[j](pt[1])*dpsi[k](pt[2])

# This is for D2
#for k in range(3):
#    for j in range(3):
#        for i in range(3):
#            for pt in pts:
#                # (2,0,0)
#                print ddpsi[i](pt[0])*psi[j](pt[1])*psi[k](pt[2])
#                # (1,1,0)
#                print dpsi[i](pt[0])*dpsi[j](pt[1])*psi[k](pt[2])
#                # (1,0,1)
#                print dpsi[i](pt[0])*psi[j](pt[1])*dpsi[k](pt[2])
#                # (0,2,0)
#                print psi[i](pt[0])*ddpsi[j](pt[1])*psi[k](pt[2])
#                # (0,1,1)
#                print psi[i](pt[0])*dpsi[j](pt[1])*dpsi[k](pt[2])
#                # (0,0,2)
#                print psi[i](pt[0])*psi[j](pt[1])*ddpsi[k](pt[2])

# This is for D3
#for k in range(3):
#    for j in range(3):
#        for i in range(3):
#            for pt in pts:
#                # (3,0,0)
#                print 0.0
#                # (2,1,0)
#                print ddpsi[i](pt[0])*dpsi[j](pt[1])*psi[k](pt[2])
#                # (2,0,1)
#                print ddpsi[i](pt[0])*psi[j](pt[1])*dpsi[k](pt[2])
#                # (1,2,0)
#                print dpsi[i](pt[0])*ddpsi[j](pt[1])*psi[k](pt[2])
#                # (1,1,1)
#                print dpsi[i](pt[0])*dpsi[j](pt[1])*dpsi[k](pt[2])
#                # (1,0,2)
#                print dpsi[i](pt[0])*psi[j](pt[1])*ddpsi[k](pt[2])
#                # (0,3,0)
#                print 0.0
#                # (0,2,1)
#                print psi[i](pt[0])*ddpsi[j](pt[1])*dpsi[k](pt[2])
#                # (0,1,2)
#                print psi[i](pt[0])*dpsi[j](pt[1])*ddpsi[k](pt[2])
#                # (0,0,3)
#                print 0.0
                
# This is for D4                                
for k in range(3):
    for j in range(3):
        for i in range(3):
            for pt in pts:
                # (4,0,0)
                print 0.0
                # (3,1,0)
                print 0.0
                # (3,0,1)
                print 0.0
                # (2,2,0)
                print ddpsi[i](pt[0])*ddpsi[j](pt[1])*psi[k](pt[2])
                # (2,1,1)
                print ddpsi[i](pt[0])*dpsi[j](pt[1])*dpsi[k](pt[2])
                # (2,0,2)
                print ddpsi[i](pt[0])*psi[j](pt[1])*ddpsi[k](pt[2])
                # (1,3,0)
                print 0.0
                # (1,2,1)
                print dpsi[i](pt[0])*ddpsi[j](pt[1])*dpsi[k](pt[2])
                # (1,1,2)
                print dpsi[i](pt[0])*dpsi[j](pt[1])*ddpsi[k](pt[2])
                # (1,0,3)
                print 0.0
                # (0,4,0)
                print 0.0
                # (0,3,1)
                print 0.0
                # (0,2,2)
                print psi[i](pt[0])*ddpsi[j](pt[1])*ddpsi[k](pt[2])
                # (0,1,3)
                print 0.0
                # (0,0,4)
                print 0.0
                
                
