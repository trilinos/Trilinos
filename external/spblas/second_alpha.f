      subroutine second_alpha(t_cpu)
      real*4 etime, tarray(2)
      double precision t_cpu
c
      t_cpu = etime(tarray)
      t_cpu = tarray(1) + tarray(2)
      return
      end
