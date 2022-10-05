C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      subroutine wheneq(n,array,inc,target,index,nval)
C see WHENEQ (3SCI) Math and Scientific Library, Cray Research, Inc
      DIMENSION ARRAY(*), INDEX(*)
      INA = 1
      NVAL = 0
      IF(INC .LT. 0) INA = (-INC) * (N-1) + 1
      DO 100 I=1,N
          IF( ARRAY(INA) .EQ. TARGET)THEN
            NVAL = NVAL + 1
            INDEX(NVAL) = I
          ENDIF
          INA = INA + INC
 100  CONTINUE
      RETURN
      END

