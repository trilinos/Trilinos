C-----------------------------------------------------------------------
C           --------<   precond.h   >--------
C-----------------------------------------------------------------------
c       double precision   ID(0:m2max+buff, 0:0)
        double precision   B1(0:m2max+buff,-1:1)
c       double precision B1R1(0:m2max+buff,-2:2)
c       double precision B1R2(0:m2max+buff,-3:3)
c       double precision B1R3(0:m2max+buff,-4:4)
c       double precision B1R4(0:m2max+buff,-5:5)
c       double precision B1R5(0:m2max+buff,-6:6)
c       double precision   B2(0:m2max+buff,-2:2)
        double precision B2R1(0:m2max+buff,-3:3)
c       double precision B2R2(0:m2max+buff,-4:4)
c       double precision B2R3(0:m2max+buff,-5:5)
c       double precision B2R4(0:m2max+buff,-6:6)
c       double precision B2R6(0:m2max+buff,-8:8)
        double precision   R1(0:m2max+buff,-1:1)
c       double precision   R2(0:m2max+buff,-2:2)
c       double precision   R3(0:m2max+buff,-3:3)
c       double precision   R4(0:m2max+buff,-4:4)
c       double precision   R5(0:m2max+buff,-5:5)
c       double precision   R6(0:m2max+buff,-6:6)
C-----------------------------------------------------------------------
        common/prec1/ R1, B1, B2R1
c       common/prec1/B1,B1R1,B1R2,B1R3,B1R4,B1R5,B2,B2R1,B2R2,B2R3,B2R4
c    _              ,R1,R2,R3,R4,R5,R6,ID
C-----------------------------------------------------------------------
