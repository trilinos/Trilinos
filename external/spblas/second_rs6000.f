      subroutine second_rs6000(t_cpu)
      integer mclock,itimer
      double precision t_cpu
c
      itimer = mclock()
      t_cpu = itimer
      t_cpu = t_cpu/100.0
      return
      end
