C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INIGRD (FROM, TO, GRAD, NINT, NPTS, ARRAY)
C=======================================================================

C     INIGRD: Initialize array ARRAY with NPTS values.
C             values are calculated to have a gradient of GRAD
C             ranging from FROM to TO with NINT segments.

C --- FROM - IN - Minimum value of range
C --- TO   - IN - Maximum value of range
C --- GRAD - IN - Gradient
C --- NINT - IN - Number of intervals in range
C --- NPTS - IN - Number of points in range.  Can be less than NINT if
C                 do not need full range, but want gradient spacing
C                 based on full range.  Normally NPTS = NINT + 1
C --- ARRAY- OUT- Range of values

      REAL    ARRAY(NPTS)
      LOGICAL NOGRAD

      IF (FROM .EQ. TO  .OR. NINT .LE. 0) THEN
         CALL PRTERR ('PROGRAM',
     *      'invalid values passed to INIGRD')
         RETURN
      END IF

      NOGRAD = (ABS(GRAD - 1.0)/NINT .LE. 1.0e-7)

      IF (NOGRAD) THEN
         D3 = 1.0 / NINT
         DO 10 I=1, NPTS
            ARRAY(I) = D3 * (I-1)
   10    CONTINUE
      ELSE
         ARRAY(1) = 0.0
         D3 = (1.0 - GRAD) / (1.0 - GRAD**NINT)
         DO 20 I=2,NPTS
            ARRAY(I) = ARRAY(I-1) + D3 * GRAD**(I-2)
   20    CONTINUE
      END IF

      DO 30 I=1, NPTS
         ARRAY(I) = FROM + ARRAY(I) * (TO - FROM)
   30 CONTINUE

      RETURN
      END
