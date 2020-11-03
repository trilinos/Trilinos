C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE QALAB (DBORD0, DVIEW0, CHLSIZ, DOQA, DOAXIS, DOLAB,
     &   CAPTN, TITLE, CREATE, MODIFY, DRAW, DLEGND, BLKCOL,  *)
C=======================================================================

C   --*** QALAB *** (BLOT) Draw standard label (PLT)
C   --   Written by Amy Gilkey - revised 02/02/88
C   --
C   --QALAB starts a "standard" plot by:
C   --   o starting a new plot,
C   --   o outlining the display area (if DOQA true),
C   --   o displaying the database title and creation, modification and
C   --        plot information (if DOQA true),
C   --   o displaying the plot caption, and
C   --   o returning the legend coordinates.
C   --
C   --Parameters:
C   --   DBORD0 - IN - the plot boundary, including legend
C   --      (left, right, bottom, top)
C   --   DVIEW0 - IN - the plot view boundary (left, right, bottom, top)
C   --   CHLSIZ - IN - the size of a character line
C   --   DOQA - IN - true iff QA information is to be included in label
C   --   DOAXIS - IN - true iff axis is to be numbered (not in this routine)
C   --   DOLAB - IN - true iff axis is to be labeled (not in this routine)
C   --   CAPTN - IN - the three-line plot caption
C   --   TITLE - IN - the database title
C   --   CREATE - IN - the database creation code name, version, date, time
C   --   MODIFY - IN - the database modification code name, version, date, time
C   --   DRAW - IN - the database plot code name, version, date, time
C   --   DLEGND - OUT - the location of the legend (device units)
C   --      (left, right, bottom, top)
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - the return statement if the cancel function is active

C   --Routines Called:
C   --   GRABRT - (GRPLIB) Check for plot set abort
C   --   GRBOX  - (GRPLIB) Draw box
C   --   GRPBEG - (PLTLIB) Begin a new plot
C   --   GRTEXT - (GRPLIB) Display a software/hardware string
C   --   GRYCEN - (GRPLIB) Find center lines of text area
C   --   LENSTR - (STRLIB) Find string length

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      include 'params.blk'
      include 'dbnums.blk'

      REAL DBORD0(KTOP), DVIEW0(KTOP)
      REAL CHLSIZ
      LOGICAL DOQA, DOAXIS, DOLAB
      CHARACTER*80 CAPTN(3)
      CHARACTER*80 TITLE
      CHARACTER*(MXSTLN) CREATE(4), MODIFY(4), DRAW(4)
      REAL DLEGND(KTOP)
      CHARACTER*256 CDUM1(5)
      CHARACTER*8 CDUM2(10)
      INTEGER IDUM2, IDUM3(10)
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL GRABRT

C   --Set background color

      CALL SETBCK (2, CDUM1, IDUM2, IDUM3, CDUM2, *100)
      GOTO 110
  100 CONTINUE
      WRITE (*, *) 'Problem setting background color'

  110 CONTINUE

C   --Begin a new plot

      CALL GRPBEG

C   --Set foreground color

      CALL UGRCOL (0, BLKCOL)

C   --Set up layout

      DLEGND(KLFT) = DVIEW0(KRGT) + CHLSIZ
      DLEGND(KRGT) = DBORD0(KRGT)
      DLEGND(KTOP) = DVIEW0(KTOP)
      DLEGND(KBOT) = DVIEW0(KBOT)

      IF (DOQA) THEN

C      --Outline the display area

         IF (GRABRT ()) RETURN 1
         CALL GRBOX ('L',
     &      DBORD0(KLFT), DBORD0(KRGT), DBORD0(KBOT), DBORD0(KTOP))

C *** Database title and dates ***

C      --Title

         IF (GRABRT ()) RETURN 1
         DXTITL = 0.5 * (DBORD0(KRGT) - DBORD0(KLFT))
         DYTITL = DBORD0(KTOP) -
     &      0.5 * (DBORD0(KTOP) - DVIEW0(KTOP) + CHLSIZ)

C      --Left-justify title
         DO 120 N1 = 1, 79
            IF (TITLE(N1:N1) .NE. ' ') GO TO 130
  120    CONTINUE
  130    CONTINUE

         CALL GRTEXC (DXTITL, DYTITL, TITLE(N1:))

C      --Creator
         call qaleg(dlegnd, chlsiz, create, 1.0, 'Cre')

C      --Modifier
         call qaleg(dlegnd, chlsiz, modify, 1.5, 'Mod')

C      --User
         call qaleg(dlegnd, chlsiz, draw,   1.5, 'Drw')

      end if
C *** Plot caption ***

      DO 140 IEND = 3, 1, -1
         IF (CAPTN(IEND) .NE. ' ') GOTO 150
  140 CONTINUE
  150 CONTINUE
      IF (IEND .GT. 0) THEN
         DXCAPT = DVIEW0(KLFT) + 0.5 * (DVIEW0(KRGT) - DVIEW0(KLFT))
         DTOP = DVIEW0(KBOT) - CHLSIZ
         IF (DOAXIS) THEN
            DTOP = DTOP - 3.0*CHLSIZ
         ELSE IF (DOLAB) THEN
            DTOP = DTOP - 1.5*CHLSIZ
         END IF
         DBOT = DBORD0(KBOT)
         CALL GRYCEN (CHLSIZ, DTOP, DBOT, IEND, IDUM)
         DO 160 I = 1, IEND
            IF (GRABRT ()) RETURN 1
            CALL GRTEXC (DXCAPT, DTOP, CAPTN(I))
            DTOP = DTOP - CHLSIZ
  160    CONTINUE
      END IF

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN
      END

      subroutine qaleg(dlegnd, chlsiz, entry, factor, label)
      include 'params.blk'
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)
      REAL DLEGND(KTOP)
      REAL CHLSIZ
      CHARACTER*(MXSTLN) entry(4)
      character*3 label

      DLEGND(KTOP) = DLEGND(KTOP) - factor * CHLSIZ
C ... Limit creator, modified by, and drawn by strings to 20 characters
C     Larger strings cause problems with metafiles
      L1 = MIN(20, LENSTR(ENTRY(1)))
      CALL GRTEXT (DLEGND(KLFT), DLEGND(KTOP),
     &     label//': '//ENTRY(1)(:L1))
      DLEGND(KTOP) = DLEGND(KTOP) - CHLSIZ

      l3 = lenstr(entry(3))
      l4 = lenstr(entry(4))
      CALL GRTEXT (DLEGND(KLFT), DLEGND(KTOP),
     &     '  '//ENTRY(3)(:l3) //
     *     '  '//ENTRY(4)(:l4))

      return
      end

