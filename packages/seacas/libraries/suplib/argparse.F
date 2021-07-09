C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      integer function argument_count()
      integer*4 arg_count4
      arg_count4 = iargc()
      argument_count = arg_count4
      return
      end

      subroutine get_argument(which, optvalue, length)
      integer which
      character*(*) optvalue
      integer length
#if defined(__XLF__) || defined(__INTEL_COMPILER)
      integer*4 which4
      which4 = which
      call getarg(which4, optvalue)
#else
      call getarg(which, optvalue)
#endif
      length = lenstr(optvalue)
      return
      end
