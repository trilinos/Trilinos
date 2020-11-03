C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SRTSBC (MXSFLG, NPSBC, NPELEM, NNXK, NXK, NSFLG, NSLEN,
     &   NSPTR, NVLEN, NVPTR, LSTSBC, NELEMS, NSIDEN, NSIDES, NNSBC,
     &   NSLIST, NVLIST, NBCSID)
C***********************************************************************

C  SUBROUTINE SRTSBC = SORTS THE LIST OF SIDE BOUNDARY CARDS

C***********************************************************************

C  VARIABLES USED:
C     NSFLG  = THE ARRAY OF FLAG VALUES
C     NSLEN  = NUMBER OF ELEMENTS IN ILIST ASSOCIATED WITH EACH FLAG
C     NSPTR  = POINTER TO THE FIRST ELEMENT IN LIST FOR EACH FLAG
C     ILIST  = THE ELEMENT LIST
C     KKK    = THE NUMBER OF ELEMENTS IN THE MESH
C     MXSFLG = THE NUMBER OF ENTRIES IN THE BOUNDARY LIST
C     FOUND  = .TRUE. IF A NEW UNIQUE FLAG HAS BEEN FOUND

C***********************************************************************

      DIMENSION NXK (NNXK, NPELEM)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG), NSPTR (MXSFLG)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG)
      DIMENSION LSTSBC (NPSBC), NELEMS (NPSBC)
      DIMENSION NSIDES (NPSBC), NSIDEN (NPSBC)

      LOGICAL FOUND

      IFLAG  = -1
      NSLIST = 0
      NBCSID = 0
      IBEGIN = 1

  100 CONTINUE
      FOUND = .FALSE.

      DO 110 I = IBEGIN, NNSBC, 3
         IF (LSTSBC (I) .LT. 0) THEN
            IF (FOUND) THEN
               IF (IFLAG .EQ. ABS (LSTSBC (I))) THEN
                  NSLIST = NSLIST + 1
                  NELEMS (NSLIST) = LSTSBC (I + 1)
                  NSIDES (NSLIST) = LSTSBC (I + 2)
                  NSLEN (NBCSID) = NSLEN (NBCSID) + 1
                  LSTSBC (I) = 0
               ENDIF
            ELSE
               FOUND = .TRUE.
               NBCSID = NBCSID + 1
               IFLAG =  - LSTSBC (I)
               NSFLG (NBCSID) = IFLAG
               NSLEN (NBCSID) = 1
               NSLIST = NSLIST + 1
               NSPTR (NBCSID) = NSLIST
               NELEMS (NSLIST) = LSTSBC (I + 1)
               NSIDES (NSLIST) = LSTSBC (I + 2)
               LSTSBC (I) = 0
               IBEGIN = I
            ENDIF
         ENDIF
  110 CONTINUE

      IF (FOUND) THEN
         GOTO 100
      ELSE

C  PUT ALL THE NODES ATTACHED TO THE ELEMENT BCC INTO THE
C  NSIDEN LIST

         NVLIST = 0
         DO 130 I = 1, NBCSID
            ISTART = NSPTR (I)
            IEND = NSPTR (I) + NSLEN (I) - 1
            NVPTR (I) = NVLIST + 1
            DO 120 J = ISTART, IEND
               J1 = NSIDES (J)
               J2 = J1 + 1
               IF (J2 .EQ. 5) J2 = 1
               NVLIST = NVLIST + 1
               NSIDEN (NVLIST) = NXK (J1, NELEMS (J))
               NVLIST = NVLIST + 1
               NSIDEN (NVLIST) = NXK (J2, NELEMS (J))
  120       CONTINUE
            NVLEN (I) = NVLIST - NVPTR (I) + 1
  130    CONTINUE
         RETURN
      ENDIF

      END
