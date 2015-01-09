# to use this, do 'module load percept' first
import time

import sys

from math import *
from random import *
import numpy
import subprocess
import time
import matplotlib.pyplot as plt

#########################################################################################################################
tt='1:07'
tt='07.810'

def convert_time(tt):
    tt0=tt.split(':')
    ttsec=0.0
    if tt0.__len__() == 3:
        ttsec=float(tt0[0])*60.0*60.0+float(tt0[1])*60.0+float(tt0[2])
    elif tt0.__len__() == 2:
        ttsec=float(tt0[0])*60.0+float(tt0[1])
    else:
        ttsec=float(tt0[0])

    return ttsec


timestest=['1:07', '7', '7.2', '1:07.2', '02:01:07.2']
for tt in timestest:
    print 'test ', tt, " = ", convert_time(tt)

#########################################################################################################################

# runs conchas for 1,2,4,8 procs for a strong scaling study

#feti = True
#fetistr = 'feti'
nRR = 'R'
doRun = False
#doRun = True

nprocs=[1,2,4,8]
#nprocs=[4,8]
rtime=[]
ii = 0
timingstrings=['Perf: RunSierra::Domain::execute', 'Matrix Assemble', 'Matrix Solve', 'Perf_ppe', 'Perf_momentum', 'Perf_ppe_solve', 'Perf_mom_solve', 'Perf_ppe_asmb', 'Perf_mom_asmb',
               'Perf_bicg_solve', 'Perf_mom_apply', 'Load Complete', 'Perf_bicg_load_complete']
timingstrings=['Perf_ppe', 'Perf_momentum', 'Perf_ppe_solve', 'Perf_mom_solve',
               'Perf_ppe_asmb', 'Perf_mom_asmb',
               'Perf_ppe_asmb_cont', 'Perf_mom_asmb_cont',
               'Perf_ppe_asmb_bc', 'Perf_mom_asmb_bc',
               'Perf_ppe_asmb_comp', 'Perf_mom_asmb_comp']

timingstrings=['Perf_bicg_load_complete',
               'Perf_bicg_resNorm',
'Perf_norm',
'Perf_norm_1',
'Perf_norm_2',
'Perf_bicg_par_misc',
'Perf_bicg_par_mvprod',
'Perf_dot',
'Perf_dot_1',
'Perf_dot_2',
'Perf_bicg_precond',
'Perf_ParAssemb::run',
'Perf_ParAssemb::pack',
'Perf_ParAssemb::allocateBuffers',
'Perf_ParAssemb::comm',
'Perf_ParAssemb::unpack_and_resend',
'Perf_ParAssemb::u_and_r_unpack',
'Perf_ParAssemb::u_and_r_resend',
'Perf_ParAssemb::u_and_r_unpack2']

timingstrings=['Perf_bicg_load_complete',
               'Perf_bicg_load_complete_0',
               'Perf_bicg_load_complete_1',
               'Perf_bicg_load_complete_2',
               'Perf_bicg_load_complete_3',
               'Perf_bicg_load_complete_4',
               'Perf_bicg_load_complete_5',
               'Perf_bicg_load_complete_6']

timingstrings=['Perf_bicg_load_complete',
               'Perf_bicg_resNorm',
'Perf_norm',
'Perf_norm_1',
'Perf_norm_2',
'Perf_bicg_par_misc',
'Perf_bicg_par_mvprod',
'Perf_dot',
'Perf_dot_1',
'Perf_dot_2',
'Perf_bicg_precond',
'Perf_ParAssemb::run']

timingstrings=['Perf_ppe', 'Perf_momentum', 'Perf_ppe_solve', 'Perf_mom_solve', 'Perf_ppe_asmb', 'Perf_mom_asmb']

timingstrings=['Perf: RunSierra::Domain::execute', 'Matrix Assemble', 'Matrix Solve', 'Perf_ppe', 'Perf_momentum', 'Perf_ppe_solve', 'Perf_mom_solve', 'Perf_ppe_asmb', 'Perf_mom_asmb', 'Load Complete']


timingstrings1=[]

for tstr in timingstrings:
    rtime.append([])
    timingstrings1.append(tstr.replace(' ','_'))

for nproc in nprocs:

  # copy the right grid to a proxy grid name
  #  res = subprocess.check_output("cp slantTheta%s.g slant.g" % theta, shell=True)

  # run conchas using a bash script
  print ""
  print " nproc=  \n" , nproc
  if doRun: res = subprocess.check_output("./run-par-N %d" % nproc, shell=True)
  print ""

  # save the results
  if doRun: res = subprocess.check_output("cp edgeOpenJet%s.log edgeOpenJet%s-%d.log" % (nRR, nRR, nproc), shell=True)

  jj = 0
  for tstr in timingstrings:
      res = ''
      try:
          res = subprocess.check_output("grep -s '%s' edgeOpenJet%s-%d.log" % (tstr, nRR, nproc), shell=True)
      except subprocess.CalledProcessError:
          print "grep error: couldn't find timing string: %s in edgeOpenJet%s-%d.log" % (tstr, nRR, nproc)
          res = ''
          
      tt2 = 0
      if len(res) > 0:
          res = res.replace(tstr, timingstrings1[jj])
          r_split = res.split()
          print "time = res.split()= " , r_split
          tt2 = convert_time(r_split[2])/float(nproc)
      rtime[jj].append( tt2 )
      jj += 1

  ii += 1

# if no data at 1 proc, approximate from nproc=2 value
jj = 0
for tstr in timingstrings:
    if rtime[jj][0] == 0:
        rtime[jj][0] = 2.0*rtime[jj][1]
    jj += 1

print "rtime= " , rtime

ax1 = plt.subplot(121)
#ax1 = plt.subplot()
#ax1.set_xscale("log")
#ax1.set_yscale("log")
ax1.set_xlabel("nprocs")
ax1.set_ylabel("sec")
#ax1.set_xlim(1e1, 1e3)
#ax1.set_ylim(1e2, 1e3)
#ax1.set_aspect(1)
ax1.set_title("strong scaling")

#symstrs=['o-','+-','x-','*-','^-','#-']
symstrs=['o-', 's-', 'p-', '+-', 'x-', '*-', 'v-', '^-', 'h-', 'H-', 'D-', '<-', '>-', '1-', '2-', '3-', '4-',  'd-', '|-', '_-', '.-']
lensym=len(symstrs)
for ii in range(0,len(timingstrings1)):
    print " "
    print timingstrings1[ii], " ", rtime[ii]
    ax1.plot(nprocs, rtime[ii], symstrs[(ii+1) % lensym], label=timingstrings1[ii])
ax1.legend(loc='upper left')

ax2 = plt.subplot(122)
#ax2 = plt.subplot()
#ax2.set_xscale("log")
#ax2.set_yscale("log")
ax2.set_xlabel("nprocs")
ax2.set_ylabel("speedup")
#ax2.set_xlim(1e1, 1e3)
#ax2.set_ylim(1e2, 1e3)
#ax2.set_aspect(1)
ax2.set_title("strong scaling - speedup")

nprocs1=[]
for ii in range(0,len(nprocs)):
    nprocs1.append( float(nprocs[ii])/float(nprocs[0]) )

print "nprocs= ", nprocs, " nprocs1= " , nprocs1

ax2.plot(nprocs, nprocs1, symstrs[0], label="ideal")
for ii in range(0,len(timingstrings1)):
    stime=[]
    for jj in range(0,len(nprocs)):
        rt = rtime[ii][jj]
        if rt == 0: rt = 1.e-12
        stime.append(rtime[ii][0]/rt)
    print " "
    print timingstrings1[ii], " ", stime

    ax2.plot(nprocs, stime, symstrs[(ii+1) % lensym], label=timingstrings1[ii])
ax2.legend(loc='upper left')

plt.draw()
plt.show()

#  start = time.time()
#  result = l2Norm.evaluate(ff_Tnd)
#  print "|ff_Tnd|= " , result, " time= ", time.time()-start

#  start = time.time()
#  result = l2Norm.evaluate(sf_exact)
#  print "|exact|= " , result, " time= ", time.time()-start

#test_result = l2Norm.evaluate(sf_exact)
#ff_result = l2Norm.evaluate(ff_Tnd)
#expected_test_result = 48.*sqrt(3041./5.)  # result from Mathematica

#print "norm exact= " , test_result, " expected_test_result= " , expected_test_result, " FEM norm = " , ff_result


#########################################################################################################################
