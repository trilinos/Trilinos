
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
      integer*4 which4

      which4 = which
      call getarg(which4, optvalue)
      length = lenstr(optvalue)
      return
      end
