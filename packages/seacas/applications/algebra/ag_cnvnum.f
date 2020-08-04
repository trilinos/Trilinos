C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CNVNUM (NUMSTR, NUM, IERR)
C=======================================================================

C   --*** CNVNUM *** (ALGEBRA) Convert an integer string to a number
C   --   Written by Amy Gilkey - revised 02/28/86
C   --
C   --CNVNUM converts an integer string into a number.
C   --
C   --Parameters:
C   --   NUMSTR - IN - the integer string (<= 10 digits)
C   --   NUM - OUT - the number
C   --   IERR - OUT - 0 iff no error occurred

      CHARACTER*(*) NUMSTR
      CHARACTER*10 INTSTR

      IERR = 0

      L = INDEX (NUMSTR, ' ') - 1
      IF (L .LT. 0) L = LEN(NUMSTR)
      IF ((L .LE. 0) .OR. (L .GT. 10)) THEN
         IERR = 1
         GOTO 100
      END IF
      INTSTR = ' '
      INTSTR(10-L+1:) = NUMSTR(:L)

      READ (INTSTR, '(I10)', IOSTAT=IERR) NUM

  100 CONTINUE
      RETURN
      END
