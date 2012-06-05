C $Id: srtsbc.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: srtsbc.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:16:37  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:16:36  gdsjaar
c Initial revision
c 
C
CC* FILE: [.RENUM]SRTSBC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SRTSBC (MXSFLG, NPSBC, NPELEM, NNXK, NXK, NSFLG, NSLEN,
     &   NSPTR, NVLEN, NVPTR, LSTSBC, NELEMS, NSIDEN, NSIDES, NNSBC,
     &   NSLIST, NVLIST, NBCSID)
C***********************************************************************
C
C  SUBROUTINE SRTSBC = SORTS THE LIST OF SIDE BOUNDARY CARDS
C
C***********************************************************************
C
C  VARIABLES USED:
C     NSFLG  = THE ARRAY OF FLAG VALUES
C     NSLEN  = NUMBER OF ELEMENTS IN ILIST ASSOCIATED WITH EACH FLAG
C     NSPTR  = POINTER TO THE FIRST ELEMENT IN LIST FOR EACH FLAG
C     ILIST  = THE ELEMENT LIST
C     KKK    = THE NUMBER OF ELEMENTS IN THE MESH
C     MXSFLG = THE NUMBER OF ENTRIES IN THE BOUNDARY LIST
C     FOUND  = .TRUE. IF A NEW UNIQUE FLAG HAS BEEN FOUND
C
C***********************************************************************
C
      DIMENSION NXK (NNXK, NPELEM)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG), NSPTR (MXSFLG)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG)
      DIMENSION LSTSBC (NPSBC), NELEMS (NPSBC)
      DIMENSION NSIDES (NPSBC), NSIDEN (NPSBC)
C
      LOGICAL FOUND
C
      IFLAG  = -1
      NSLIST = 0
      NBCSID = 0
      IBEGIN = 1
C
  100 CONTINUE
      FOUND = .FALSE.
C
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
C
      IF (FOUND) THEN
         GOTO 100
      ELSE
C
C  PUT ALL THE NODES ATTACHED TO THE ELEMENT BCC INTO THE
C  NSIDEN LIST
C
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
C
      END
