C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CKNONE (NVAL, ISSEL, VALNAM, *)
C=======================================================================

C   --*** CKNONE *** (ETCLIB) Check number of values is zero
C   --   Written by Amy Gilkey - revised 12/23/87
C   --
C   --CKNONE prints an error message if the number of values is zero.
C   --
C   --Parameters:
C   --   NVAL - IN - the value being checked
C   --   ISSEL - IN - print none selected error message iff true
C   --   VALNAM - IN - the name of the value being checked (plural)
C   --   * - return statement if error

      LOGICAL ISSEL
      CHARACTER*(*) VALNAM

      CHARACTER*1024 ERRMSG

      IF (NVAL .LE. 0) THEN
         IF (ISSEL) THEN
            ERRMSG = 'No ' // VALNAM // ' are selected'
            CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
         ELSE
            ERRMSG = 'There are no ' // VALNAM
            CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
         END IF
         RETURN 1
      END IF

      RETURN
      END
