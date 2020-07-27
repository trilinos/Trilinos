C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LTSORT (MDIM, LINK, ID, IPNTR, ADDLNK)
C***********************************************************************

C  SUBROUTINE LTSORT = LOOKUP TABLE SORT FOR DATA POINTER ARRAYS

C***********************************************************************

C  VARIABLES USED:
C     MDIM   = DIMENSION OF LINK ARRAY,  AND BASE FOR LOOKUP START
C     LINK   = LOOKUP TABLE ARRAY OF ID'S AND POINTERS
C              LINK (1, I) = ID VALUE STORED IN I'TH ROW  (0 IF EMPTY)
C              LINK (2, I) = DATA POINTER ASSOCIATED W/THIS I'TH ID VALUE
C     ID     = THE ID OF THE DATA BEING FOUND OR PLACED
C     IPNTR  = THE DATA POINTER ASSOCIATED WITH THE ID BEING USED
C     ADDLNK = .TRUE. IF DATA IS BEING PLACED IN THE LOOKUP TABLE
C            = .FALSE. IF DATA IS BEING FOUND ONLY

C***********************************************************************

      DIMENSION LINK (2, MDIM)

      LOGICAL ADDLNK

C  CALCULATE THE BEGINNING LOOKUP VALUE

      if (id .lt. 0) stop 'LTSORT: Internal error'

      HOLD = DBLE(ID) * 3.1830989
      LOOKUP =  INT((HOLD - INT (HOLD)) * DBLE(MDIM) + 1)

C  SET UP THE LOOP TO ONLY SEARCH THROUGH THE TABLE ONCE

      DO 100 I = 1, MDIM

C  IF LOOKUP SPOT IS EMPTY THEN FILL AND RETURN IF ADDING AND IPNTR .NE. 0
C  OR FLAG IPNTR AS BEING EMPTY AND RETURN IF FINDING

         IF (LINK (1, LOOKUP) .EQ. 0) THEN
            IF ( (ADDLNK) .AND. (IPNTR .NE. 0)) THEN
               LINK (1, LOOKUP) = ID
               LINK (2, LOOKUP) = IPNTR
            ELSEIF (IPNTR .NE. 0) THEN
               IPNTR = 0
            ENDIF
            RETURN

C  IF LOOKUP SLOT IS FULL,  CHECK TO SEE IF IT MATCHES THE CURRENT ID
C  IF IT MATCHES AND IF ADDING,  SET THE NEW POINTER  (OVERWRITE)
C  IF IT MATCHES AND IF FINDING,  RETURN THE CORRECT POINTER
C  IF NO MATCH,  THEN INCREMENT LOOKUP AND TRY AGAIN IN THE TABLE

         ELSE
            IF (ID .EQ. LINK (1, LOOKUP)) THEN
               IF (ADDLNK) THEN
                  LINK (2, LOOKUP) = IPNTR
               ELSE
                  IPNTR = LINK (2, LOOKUP)
               ENDIF
               RETURN
            ELSE
               LOOKUP = LOOKUP + 1
               IF (LOOKUP .GT. MDIM)LOOKUP = 1
            ENDIF
         ENDIF
  100 CONTINUE

C  ACT ON THE EXHAUSTED SEARCH

      IF (ADDLNK) THEN
         CALL MESAGE ('LOOKUP TABLE OVERFLOW')
         CALL MESAGE ('SERIOUS DATA PROBLEMS HAVE BEEN CAUSED')
      ELSE
         IPNTR = 0
      ENDIF
      RETURN

      END
