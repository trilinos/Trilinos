      double precision function second()
      integer mclock,itimer
c
      itimer = mclock()
      second = itimer
      second = second/100.0D0
      return
      end
