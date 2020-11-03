C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LOGERR (ERRTYP, ERRMSG, IUN)
C=======================================================================
C   --*** LOGERR *** (ETCLIB) Print error message to output file
C   --
C   --Parameters:
C   --   ERRTYP - IN - the type of error:
C   --      'FATAL', 'PROGRAM', 'ERROR', 'WARNING',
C   --      'CMDERR', 'CMDWARN', 'CMDREQ', 'CMDSPEC'
C   --   ERRMSG - IN - the error message
C   --   IUN - IN - the output unit number

      CHARACTER*(*) ERRTYP
      CHARACTER*(*) ERRMSG

      IF (ERRTYP .EQ. 'FATAL') THEN
         WRITE (IUN, *)
         WRITE (IUN, 10) ERRMSG
   10    FORMAT (' FATAL ERROR - ', A)
      ELSE IF (ERRTYP .EQ. 'PROGRAM') THEN
         WRITE (IUN, *)
         WRITE (IUN, 20) ERRMSG, ' - email code sponsor'
   20    FORMAT (' PROGRAM ERROR - ', A, A)
      ELSE IF (ERRTYP .EQ. 'ERROR') THEN
         WRITE (IUN, *)
         WRITE (IUN, 30) ERRMSG
   30    FORMAT (' ERROR - ', A)
      ELSE IF (ERRTYP .EQ. 'WARNING') THEN
         WRITE (IUN, *)
         WRITE (IUN, 40) ERRMSG
   40    FORMAT (' WARNING - ', A)
      ELSE IF (ERRTYP .EQ. 'CMDERR') THEN
         WRITE (IUN, 50) ERRMSG
   50    FORMAT (' *** ERROR - ', A)
      ELSE IF (ERRTYP .EQ. 'CMDWARN') THEN
         WRITE (IUN, 60) ERRMSG
   60    FORMAT (' *** WARNING - ', A)
      ELSE IF (ERRTYP .EQ. 'CMDREQ') THEN
         WRITE (IUN, 70) ERRMSG
   70    FORMAT (' *** ', A)
      ELSE IF (ERRTYP .EQ. 'CMDSPEC') THEN
         WRITE (IUN, 80) ERRMSG
   80    FORMAT (' *** ', A)
      END IF

      RETURN
      END
