! Command line argument functions for NAGWare f90 2.2

      integer function mpir_iargc()
      use f90_unix
      mpir_iargc = iargc()
      return
      end

      subroutine mpir_getarg( i, s )
      use f90_unix
      integer       i
      character*(*) s
      call getarg(i,s)
      return
      end
