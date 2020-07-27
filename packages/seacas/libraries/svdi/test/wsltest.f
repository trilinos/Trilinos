C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      PROGRAM WTEST05
      REAL X(102), Y(102)
      CHARACTER COLOR(7)*4

      DATA COLOR/'RED','GREE','YELL','BLUE','MAGE','CYAN','WHIT'/

      CALL WSTART(0.,0)
      CALL WTTYPE('SOFT')
      DO 10 I=1,102
         X(I) = I/100.*2.*3.141592
         Y(I) = SIN(X(I))
   10 CONTINUE
      CALL WXAXIS(.2,.2,.7,0.,13.,1.,'Sines of the Times')
      CALL WYAXIS(.2,.2,.5,-1.,1.,.5,' ')
      CALL WLNSTY('DOTT')
      CALL WGRID(0,1)
      CALL WLNSTY('SOLI')
      CALL WMARK('DIAM')
      DO 20 I=1,7
         CALL WSPACE('NDC')
         CALL WMSIZE(I/500.)
         CALL WSPACE('USER')
         CALL WCOLOR(COLOR(I))
         CALL WDRAW(X,Y,102,10)
         DO 15 J=1,102
            X(J) = X(J) + 1.
   15    CONTINUE
   20 CONTINUE
      CALL WCOLOR(COLOR(3))
      CALL WEND
      STOP
      END
