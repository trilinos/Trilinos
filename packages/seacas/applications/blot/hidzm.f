C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDZM (ZOOM, XZMMIN, XZMMAX, YZMMIN, YZMMAX,
     &   XINMIN, XINMAX, YINMIN, YINMAX, LENF, NLNKF, LINKF,
     &   XN, YN, HIDEF)
C=======================================================================

C   --*** HIDZM *** (MESH) Hide face outside zoom window
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --HIDZM marks a face as hidden if the maximum box around
C   --all its nodes is outside the zoom window.  Some faces
C   --that are outside the zoom window may not be flagged.
C   --
C   --Parameters:
C   --   ZOOM - IN/OUT - reset if no nodes outside zoom window
C   --   XZMMIN, XZMMAX, YZMMIN, YZMMAX - IN - the box enclosing the
C   --      zoom window
C   --   XINMIN, XINMAX, YINMIN, YINMAX - OUT - the box enclosing the
C   --      minimum needed window
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   XN, YN - IN - the nodal coordinates
C   --   HIDEF(i) - IN/OUT - set true iff face i is hidden
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      PARAMETER (KFVIS=0, KFNODH=10, KFPOUT=20, KFOUT=90, KFAWAY=100)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL ZOOM
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      REAL XN(*), YN(*)
      INTEGER HIDEF(*)

      ZOOM = .FALSE.

      XXMIN = XZMMIN
      XXMAX = XZMMAX
      YYMIN = YZMMIN
      YYMAX = YZMMAX

      DO 120 IELB = 1, NELBLK
         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
         DO 110 IFAC = LENF(IELB-1)+1, LENF(IELB)

C         --Calculate enclosing X-Y-Z box for face

            XMIN = XN(LINKF(IXL0+1))
            XMAX = XMIN
            YMIN = YN(LINKF(IXL0+1))
            YMAX = YMIN
            DO 100 ILINK = 2, NLNKF(IELB)
               X = XN(LINKF(IXL0+ILINK))
               Y = YN(LINKF(IXL0+ILINK))
               XMIN = MIN (XMIN, X)
               XMAX = MAX (XMAX, X)
               YMIN = MIN (YMIN, Y)
               YMAX = MAX (YMAX, Y)
  100       CONTINUE

C         --Hide face if outside zoom window

            IF ((XMAX .LE. XZMMIN) .OR. (XMIN .GE. XZMMAX)
     &         .OR. (YMAX .LE. YZMMIN) .OR. (YMIN .GE. YZMMAX)) THEN
               HIDEF(IFAC) = KFPOUT
               ZOOM = .TRUE.
            ELSE IF ((XMIN .LE. XZMMIN) .OR. (XMAX .GE. XZMMAX)
     &         .OR. (YMIN .LE. YZMMIN) .OR. (YMAX .GE. YZMMAX)) THEN
               XXMIN = MIN (XXMIN, XMIN)
               XXMAX = MAX (XXMAX, XMAX)
               YYMIN = MIN (YYMIN, YMIN)
               YYMAX = MAX (YYMAX, YMAX)
               ZOOM = .TRUE.
            END IF

            IXL0 = IXL0 + NLNKF(IELB)
  110    CONTINUE
  120 CONTINUE

      IF (.NOT. ZOOM) THEN
         XINMIN = XZMMIN
         XINMAX = XZMMAX
         YINMIN = YZMMIN
         YINMAX = YZMMAX
      ELSE
         XINMIN = XXMIN
         XINMAX = XXMAX
         YINMIN = YYMIN
         YINMAX = YYMAX
      END IF

      IF (.NOT. ZOOM) RETURN

      nout = 0
      nhid = 0
      DO 150 IELB = 1, NELBLK
         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
         DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (HIDEF(IFAC) .EQ. KFPOUT) THEN

C            --Calculate enclosing X-Y-Z box for face

               XMIN = XN(LINKF(IXL0+1))
               XMAX = XMIN
               YMIN = YN(LINKF(IXL0+1))
               YMAX = YMIN
               DO 130 ILINK = 2, NLNKF(IELB)
                  X = XN(LINKF(IXL0+ILINK))
                  Y = YN(LINKF(IXL0+ILINK))
                  XMIN = MIN (XMIN, X)
                  XMAX = MAX (XMAX, X)
                  YMIN = MIN (YMIN, Y)
                  YMAX = MAX (YMAX, Y)
  130          CONTINUE

C            --Hide face if outside zoom window

               IF ((XMAX .LE. XXMIN) .OR. (XMIN .GE. XXMAX)
     &            .OR. (YMAX .LE. YYMIN) .OR. (YMIN .GE. YYMAX)) THEN
                  HIDEF(IFAC) = KFOUT
                  nhid = nhid + 1
               ELSE
                  XINMIN = MIN (XINMIN, XMIN)
                  XINMAX = MAX (XINMAX, XMAX)
                  YINMIN = MIN (YINMIN, YMIN)
                  YINMAX = MAX (YINMAX, YMAX)
                  nout = nout + 1
               END IF
            END IF

            IXL0 = IXL0 + NLNKF(IELB)
  140    CONTINUE
  150 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1)) then
         write (*, '(1x,a,i5)') 'faces outside zoom window =', nout
         write (*, '(1x,a,i5)') 'faces outside zoom slack =', nhid
      end if

      RETURN
      END
