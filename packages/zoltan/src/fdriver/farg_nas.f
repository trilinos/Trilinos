! Command line argument functions for NASoftware FortranPlus 2.0

      integer function mpir_iargc()
      use nas_system
      mpir_iargc = iargc()
      return
      end

      subroutine mpir_getarg( i, s )
      use nas_system
      integer       i
      character*(*) s
      call getarg(i,s)
      return
      end
