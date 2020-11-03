C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDWHE (VERB, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   XN, YN, ZN, XNN, YNN, ZNN, HIDENP, XE, YE, ZE, MAPEL, MAPND, *)
C=======================================================================

C   --*** CMDWHE *** (MESH) Process locate commands
C   --   Written by Amy Gilkey - revised 02/05/88
C   --
C   --Parameters:
C   --   VERB - IN/OUT - the verbs for the SHOW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   XN, YN, ZN - IN - the coordinates of the original nodes
C   --   XNN, YNN, ZNN - IN - the coordinates of the displayed nodes
C   --   HIDENP(i) - IN - true iff node i is visible,
C   --      displayed mesh coordinates are defined if not HIDENP(i)
C   --   XE, YE, ZE - IN - the coordinates of the element centroids
C   --
C   --Common Variables:
C   --   Uses NUMNP of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'rotopt.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      REAL XN(NUMNP), YN(NUMNP), ZN(NUMNP)
      REAL XNN(NUMNP), YNN(NUMNP), ZNN(NUMNP)
      LOGICAL HIDENP(NUMNP)
      REAL XE(NUMEL), YE(NUMEL), ZE(NUMEL)
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*(MXSTLN) WORD
      LOGICAL ISON
      REAL RNUM(3)
      CHARACTER*20 RSTR(3)
      LOGICAL WINDOW
      LOGICAL FFMATC, MATSTR
      LOGICAL LDUM(1)

      IF ((VERB .EQ. 'WHAT') .OR. (VERB .EQ. 'WHAT3')) THEN
         INLINE = ' '
         WINDOW = (VERB .EQ. 'WHAT')
         VERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

         IF (MATSTR (WORD, 'NODE', 1)) THEN
            ISON = FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)
            IF (WINDOW) THEN
               CALL PICKNP ('object location', ISON, WINDOW,
     &            NDIM, NUMNP, XNN, YNN, ZNN, HIDENP,
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            INP, *100)
            ELSE
               CALL PICKNP ('object location', ISON, WINDOW,
     &            NDIM, NUMNP, XN, YN, ZN, LDUM,
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            INP, *100)
            END IF

            CALL INTSTR (1, 0, MAPND(INP), RSTR(1), L)
            WRITE (*, 10000) 'Nearest node is node ', RSTR(1)(:L)
            RNUM(1) = XN(INP)
            RNUM(2) = YN(INP)
            IF (IS3DIM) RNUM(3) = ZN(INP)
            CALL NUMSTR (3, 4, RNUM(1), RSTR, LSTR)
            WRITE (*, 10000) 'Undeformed object location is',
     &         (' ', RSTR(I)(:LSTR), I=1,3)

          ELSE IF (MATSTR (WORD, 'ELEMENT', 1)) THEN
            IF (FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)) THEN
              CALL PRTERR ('WARNING',
     *          'CURSOR option not supported for elements')
              GO TO 100
            ENDIF
            ISON = .FALSE.
            CALL PICKNP ('object location', ISON, .FALSE.,
     &        NDIM, NUMEL, XE, YE, ZE, LDUM,
     &        .TRUE., IFLD, INTYP, RFIELD,
     &        INP, *100)

            CALL INTSTR (1, 0, MAPEL(INP), RSTR(1), L)
            WRITE (*, 10000) 'Nearest object is element ', RSTR(1)(:L)
            RNUM(1) = XE(INP)
            RNUM(2) = YE(INP)
            IF (IS3DIM) RNUM(3) = ZE(INP)
            CALL NUMSTR (3, 4, RNUM(1), RSTR, LSTR)
            WRITE (*, 10000) 'Undeformed object location is',
     &         (' ', RSTR(I)(:LSTR), I=1,3)

         ELSE
            CALL PRTERR ('CMDERR', 'Expected "NODE" or "ELEMENT"')
            GOTO 100
         END IF

      ELSE IF (VERB .EQ. 'WHERE') THEN
         INLINE = ' '
         VERB = ' '
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

         IF (MATSTR (WORD, 'NODE', 1)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'node number', 0, NOD, *100)
            INP = locint(NOD, NUMNP, MAPND)
            IF ((INP .LT. 1) .OR. (INP .GT. NUMNP)) THEN
               CALL PRTERR ('CMDERR', 'Invalid node id')
               GOTO 100
            END IF

            RNUM(1) = XN(INP)
            RNUM(2) = YN(INP)
            IF (IS3DIM) RNUM(3) = ZN(INP)
            CALL NUMSTR (3, 4, RNUM(1), RSTR, LSTR)
            WRITE (*, 10000) 'Undeformed object location is',
     &         (' ', RSTR(I)(:LSTR), I=1,NDIM)
            IF (IS3DIM .AND. (.NOT. HIDENP(INP))) THEN
               RNUM(1) = XNN(INP)
               RNUM(2) = YNN(INP)
               RNUM(3) = ZNN(INP)
               CALL NUMSTR (2, 4, RNUM(1), RSTR, LSTR)
               WRITE (*, 10000) 'Displayed window location is',
     &            (' ', RSTR(I)(:LSTR), I=1,2)
            END IF

          ELSE IF (MATSTR (WORD, 'ELEMENT', 1)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'element number', 0, IEL, *100)
            INP = locint(IEL, NUMEL, MAPEL)
            IF ((INP .LT. 1) .OR. (INP .GT. NUMEL)) THEN
               CALL PRTERR ('CMDERR', 'Invalid element number')
               GOTO 100
            END IF

            RNUM(1) = XE(INP)
            RNUM(2) = YE(INP)
            IF (IS3DIM) RNUM(3) = ZE(INP)
            CALL NUMSTR (3, 4, RNUM(1), RSTR, LSTR)
            WRITE (*, 10000) 'Undeformed object location is',
     &         (' ', RSTR(I)(:LSTR), I=1,NDIM)

         ELSE
            CALL PRTERR ('CMDERR', 'Expected "NODE" or "ELEMENT"')
            GOTO 100
         END IF

      END IF

      RETURN

  100 CONTINUE
      RETURN 1
10000  FORMAT (1X, 10A)
      END
