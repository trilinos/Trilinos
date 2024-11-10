C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDEDG (HIDEF, HIDENP, LENF, NLNKE, NLNKF, LINKF,
     &   NOTSEL, IELBST, IEDSET, NEDGES, NREF, LREF, MREF,
     *  ISBACK, NAMELB)
C=======================================================================

C   --*** HIDEDG *** (MESH) Identify 3D lines on edge of visible mesh
C   --   Written by Amy Gilkey - revised 02/26/88
C   --
C   --HIDEDG finds all lines that may make up the visible edge of the mesh.
C   --"Edges" are lines of partially visible faces or lines of
C   --totally visible faces which have no matching edge.
C   --
C   --Parameters:
C   --   HIDEF(i) - IN - true iff face i is hidden
C   --   HIDENP - IN - node status (as in HIDDEN)
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   NOTSEL - IN - true iff faces with selected element blocks
C   --      are not needed
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IEDSET - OUT - the edge line set; (0) = face defining edge
C   --   NEDGES - OUT - the number of lines in the edge set
C   --   NREF - SCRATCH - the number of references to a node in the edge set;
C   --      length = NUMNPF
C   --   LREF - SCRATCH - the last edge set index of a node in the edge set;
C   --      length = NUMNPF
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses NUMNPF of /D3NUMS/

      PARAMETER (KFVIS=0, KFNODH=10, KFPOUT=20, KFOUT=90, KFAWAY=100)
      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER HIDEF(*)
      INTEGER HIDENP(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKE(NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL NOTSEL
      INTEGER IELBST(NELBLK)
      INTEGER IEDSET(0:2,*)
      INTEGER NREF(NUMNPF), LREF(NUMNPF), MREF(NUMNPF)
      LOGICAL ISBACK(NUMNPF)
      CHARACTER*(*) NAMELB(*)

      LOGICAL ANYBRI, ANYSHE

      ANYBRI = .FALSE.
      ANYSHE = .FALSE.
      DO 100 IELB = 1, NELBLK
         IF (NLNKE(IELB) .GT. 4) ANYBRI = .TRUE.
         IF (NAMELB(IELB)(:3) .EQ. 'TET') ANYBRI = .TRUE.
         IF (NLNKE(IELB) .EQ. 4 .AND. NAMELB(IELB)(:3) .NE. 'TET')
     $        ANYSHE = .TRUE.
  100 CONTINUE

      CALL INIINT (NUMNPF, 0, LREF)

      IF (.NOT. ANYBRI) GOTO 200

C   --Mark nodes in faces that point away

      CALL INILOG (NUMNPF, .FALSE., ISBACK)

      DO 130 IELB = 1, NELBLK
         IF (NLNKE(IELB) .GT. 4 .OR. NAMELB(IELB)(:3) .eq. 'TET') THEN
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .GE. KFAWAY) THEN
                  DO 110 ILINK = 1, NLNKF(IELB)
                     ISBACK(LINKF(IXL0+ILINK)) = .TRUE.
  110             CONTINUE
               END IF
               IXL0 = IXL0 + NLNKF(IELB)
  120       CONTINUE
         END IF
  130 CONTINUE

      NEDGES = 0
      nmov = 0
      DO 190 IELB = 1, NELBLK
         IF (NLNKE(IELB) .GT. 4 .OR. NAMELB(IELB)(:3) .eq. 'TET') THEN

C         --Eliminate faces of a selected element block
            IF (NOTSEL .AND. (IELBST(IELB) .GT. 0)) GOTO 190

            DO 180 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LT. KFOUT) THEN
                  IXL0 = IDBLNK (IELB, IFAC, LENF, NLNKF) - 1

C               --Eliminate faces that do not have any nodes in faces
C               --that point away

                  NBACK = 0
                  DO 140 ILINK = 1, NLNKF(IELB)
                     IF (ISBACK(LINKF(IXL0+ILINK))) NBACK = NBACK + 1
  140             CONTINUE
                  IF (NBACK .LT. 2) GOTO 180

C               --Extract face lines

                  NOLD = NEDGES
                  N2 = LINKF(IXL0+NLNKF(IELB))
                  LREFN2 = LREF(N2)
                  LSTREF = LREFN2

                  DO 170 ILINK = 1, NLNKF(IELB)
                     N1 = N2
                     N2 = LINKF(IXL0+ILINK)
                     LREFN1 = LREFN2
                     IF (ILINK .LT. NLNKF(IELB)) THEN
                        LREFN2 = LREF(N2)
                     ELSE
                        LREFN2 = LSTREF
                     END IF

                     IF (n1 .eq. n2) goto 170
                     IF (.NOT. ISBACK(N1)) GOTO 170
                     IF (.NOT. ISBACK(N2)) GOTO 170
                     IF ((HIDENP(N1) .GT. KNVIS)
     &                  .AND. (HIDENP(N2) .GT. KNVIS)) GOTO 170
                     IF (HIDENP(N1) .GE. KNHID) GOTO 170
                     IF (HIDENP(N2) .GE. KNHID) GOTO 170
                     NMIN = MIN (N1, N2)
                     NMAX = MAX (N1, N2)

C                  --Search for line in existing lines

                     DO 150 IL = MIN (LREFN1, LREFN2, NOLD), 1, -1
                       IF (IEDSET(1,IL) .EQ. NMIN .AND.
     *                     IEDSET(2,IL) .EQ. NMAX) GOTO 160
  150                CONTINUE
  160                CONTINUE

C                  --Insert new line or pointer to face if duplicate

                     IF (IL .GT. 0) THEN
                        IF (NBACK .LE. 2) IEDSET(0,IL) = IFAC
                        nmov = nmov + 1
                     ELSE
                        NEDGES = NEDGES + 1
                        IEDSET(0,NEDGES) = IFAC
                        IEDSET(1,NEDGES) = NMIN
                        IEDSET(2,NEDGES) = NMAX
                        LREF(N1) = NEDGES
                        LREF(N2) = NEDGES
                     END IF
  170             CONTINUE
               END IF
  180       CONTINUE
         END IF
  190 CONTINUE

      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'non-shell visible edge set =', nedges
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') '***removed =', nmov
  200 CONTINUE

      IF (.NOT. ANYSHE) GOTO 320

C   --Mark nodes in shell faces with any nodes hidden

      DO 210 INP = 1, NUMNPF
         IF (LREF(INP) .LE. 0) THEN
            NREF(INP) = 0
         ELSE
            NREF(INP) = 1
         END IF
  210 CONTINUE

      CALL INILOG (NUMNPF, .FALSE., ISBACK)

      DO 250 IELB = 1, NELBLK
         IF (NLNKE(IELB) .EQ. 4 .AND. NAMELB(IELB)(:3) .ne. 'TET') THEN
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 240 IFAC = LENF(IELB-1)+1, LENF(IELB)
               NHID = 0
               DO 220 ILINK = 1, NLNKF(IELB)
                  IF (HIDENP(LINKF(IXL0+ILINK)) .GT. KNVIS)
     &               NHID = NHID + 1
  220          CONTINUE
               IF (NHID .GT. 0) THEN
                  DO 230 ILINK = 1, NLNKF(IELB)
C????                     ISBACK(LINKF(IXL0+ILINK)) = .TRUE.
  230             CONTINUE
               END IF
               IXL0 = IXL0 + NLNKF(IELB)
  240       CONTINUE
         END IF
  250 CONTINUE

      NEWEDG = NEDGES + 1
      nmov = 0

      call iniint(numnpf, 0, mref)
      do i=nedges, 1, -1
        mref(iedset(1,i)) = i
        mref(iedset(2,i)) = i
      end do

      DO 310 IELB = 1, NELBLK
         IF (NLNKE(IELB) .EQ. 4 .AND. NAMELB(IELB)(:3) .ne. 'TET') THEN

C         --Eliminate faces of a selected element block
            IF (NOTSEL .AND. (IELBST(IELB) .GT. 0)) GOTO 310

            DO 300 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LT. KFOUT) THEN
                  IXL0 = IDBLNK (IELB, IFAC, LENF, NLNKF) - 1

C               --Extract face lines

                  NOLD = NEDGES
                  N2 = LINKF(IXL0+NLNKF(IELB))
                  LREFN2 = LREF(N2)
                  LSTREF = LREFN2

                  DO 290 ILINK = 1, NLNKF(IELB)
                     N1 = N2
                     N2 = LINKF(IXL0+ILINK)
                     LREFN1 = LREFN2
                     IF (ILINK .LT. NLNKF(IELB)) THEN
                        LREFN2 = LREF(N2)
                     ELSE
                        LREFN2 = LSTREF
                     END IF

                     IF ((HIDENP(N1) .GT. KNVIS)
     &                  .AND. (HIDENP(N2) .GT. KNVIS)) GOTO 290
                     IF (HIDENP(N1) .GE. KNHID) GOTO 290
                     IF (HIDENP(N2) .GE. KNHID) GOTO 290
                     NMIN = MIN (N1, N2)
                     NMAX = MAX (N1, N2)

C                  --Search for line in existing lines

                     IMIN = max(mref(n1), mref(n2),1)
                     DO 260 IL = MIN (LREFN1, LREFN2, NOLD), IMIN, -1
                       IF (IEDSET(1,IL) .EQ. NMIN .AND.
     *                     IEDSET(2,IL) .EQ. NMAX) GOTO 270
 260                 CONTINUE
 270                 CONTINUE

                     IF (IL .GT. 0) THEN
                        IF (IL .LT. NEWEDG) THEN

C                        --Leave edge matching non-shell as edge
                           CONTINUE

                        ELSE IF (ISBACK(N1) .OR. ISBACK(N2)) THEN

C                        --Adjust pointer to face if duplicate line and either
C                        --node is hidden in a shell quad
                           NHID = 0
                           DO 280 I = 1, NLNKF(IELB)
                              IF (HIDENP(LINKF(IXL0+I)) .GT. KNVIS)
     &                           NHID = NHID + 1
  280                      CONTINUE
                           IF (NHID .LE. 0) IEDSET(0,NEDGES) = IFAC
C????C should check the ZQ coordinate

                        ELSE

C                        --Delete duplicate line if both nodes visible
C                        --in all shell quads

                           IEDSET(0,IL) = IEDSET(0,NEDGES)
                           IEDSET(1,IL) = IEDSET(1,NEDGES)
                           IEDSET(2,IL) = IEDSET(2,NEDGES)
                           NREF(N1) = NREF(N1) - 1
                           NREF(N2) = NREF(N2) - 1
                           IF (NREF(N1) .LE. 0) LREF(N1) = 0
                           IF (NREF(N2) .LE. 0) LREF(N2) = 0
                           NEDGES = NEDGES - 1
                           IF (NOLD .GT. NEDGES) NOLD = NEDGES
                           nmov = nmov + 1
                        END IF
                     ELSE

C                     --Insert new line

                        NEDGES = NEDGES + 1
                        IEDSET(0,NEDGES) = IFAC
                        IEDSET(1,NEDGES) = NMIN
                        IEDSET(2,NEDGES) = NMAX
                        NREF(N1) = NREF(N1) + 1
                        NREF(N2) = NREF(N2) + 1
                        LREF(N1) = NEDGES
                        LREF(N2) = NEDGES
                        if (mref(n1) .eq. 0) mref(n1) = nedges
                        if (mref(n2) .eq. 0) mref(n2) = nedges
                     END IF
  290             CONTINUE
               END IF
  300       CONTINUE
         END IF
  310 CONTINUE

      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'shell visible edge set =',
     &   nedges - newedg + 1
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') '***removed =', nmov
  320 CONTINUE

C   --Make sure hidden node is in IEDSET(2,x)

      DO 330 IEDG = 1, NEDGES
         IF (HIDENP(IEDSET(1,IEDG)) .GT. KNVIS) THEN
            I = IEDSET(1,IEDG)
            IEDSET(1,IEDG) = IEDSET(2,IEDG)
            IEDSET(2,IEDG) = I
         END IF
  330 CONTINUE

      RETURN
      END
