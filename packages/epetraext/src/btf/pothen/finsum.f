      subroutine   finsum   ( timhrz, timesq, timvrt, hrzcmp, sqcmpn,
     $                        vrtcmp, ccmstr, rcmstr, output )

c     ==================================================================
c     ==================================================================
c     ====  finsum -- print summary of fine block triangular        ====
c     ====            decomposition                                 ====
c     ==================================================================
c     ==================================================================

c     created by john lewis, boeing computer services, sept. 17, 1990

      integer            hrzcmp, sqcmpn, vrtcmp, output

      integer            ccmstr (*), rcmstr (*)

      real               timesq, timhrz, timvrt

c     ==================================================================

      if  ( hrzcmp .gt. 0 )  then
         write (output, 60000) hrzcmp, timhrz
         call fnrsum ( 1, hrzcmp, ccmstr, rcmstr, output ) 
      endif

      if  ( sqcmpn .gt. 0 )  then
         write (output, 61000) sqcmpn, timesq
         call fnrsum ( hrzcmp + 1, hrzcmp + sqcmpn, ccmstr, rcmstr,
     $                 output ) 
      endif

      if  ( vrtcmp .gt. 0 )  then
         write (output, 62000) vrtcmp, timvrt
         call fnrsum ( hrzcmp + sqcmpn + 1, hrzcmp + sqcmpn + vrtcmp,
     $                 ccmstr, rcmstr, output ) 
      endif

      return

60000 format (/'0fine decomposition of horizontal block (hr-hc)',
     $        /'      number of connected components:', i10,
     $        /'                       time required:', 1pe10.1,
     $        /'0                    component    rows  columns' )

61000 format (/'0fine decomposition of     square block (sr-sc)',
     $        /'         number of strong components:', i10,
     $        /'                       time required:', 1pe10.1,
     $        /'0                    component    rows  columns' )

62000 format (/'0fine decomposition of   vertical block (vr-vc)',
     $        /'      number of connected components:', i10,
     $        /'                       time required:', 1pe10.1,
     $        /'0                    component    rows  columns' )

      end

