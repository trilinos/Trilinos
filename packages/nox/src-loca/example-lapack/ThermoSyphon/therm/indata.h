c      data read from input file
c
       double precision bval1,bval2, pi
       double precision prandtl, sigma, rayleigh ,  tol
       double precision ra_step, ra_final, pr_final
       integer bc_case 
       integer debug, m, n 
       character*30 oldfile, savefile
       character*1 jobvL, jobvR

c       common /indata/ kappa, tinit, tstop, deltat,bval1,bval2, pi,
c     &  prandtl, sigma, rayleigh, stol, tol,
c     &  ra_step, ra_final, pr_final,
c     &  tout, tdiag, bc_case, dparno, iparno, debug, m, n, 
c     &  order, max_order,  testmode, old_run, oldfile, savefile,
c     &  printfile, jobvL, jobvR

