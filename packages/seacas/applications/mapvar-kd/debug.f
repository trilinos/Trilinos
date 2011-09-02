      subroutine debug(routine)
      character*(*) routine
      
      include 'tapes.blk'
      include 'debg.blk'

      data timlast /0.0/
      save timlast
      
      if (idebug .le. 0) return

      call excpus(cputim)
      delta = cputim - timlast
      timlast = cputim
      write (nout, 100) routine, cputim, delta
      write (ntpout, 100) routine, cputim, delta
 100  format(5x,'Entering Routine: ',A,' at time ', 1pe10.3,
     $     ', delta = ', 1pe10.3)
      end
