
      integer function argument_count()
      argument_count = iargc()
      return
      end

      subroutine get_argument(which, optvalue, length)
      integer which
      character*(*) optvalue
      integer length

      call getarg(which, optvalue)
      length = lenstr(optvalue)
      return
      end
