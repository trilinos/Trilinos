C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MESHUP (WIDLIN, MSHLIN, MLNTYP,
     &   IELBST, LENL, LINSET, BLKCOL, XN, YN, ZN,
     &   IDELB, *)
C=======================================================================

C   --*** MESHUP *** (MESH) Plot mesh lines
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 03/29/88
C   --   D. P. Flanagan, 07/02/82
C   --
C   --MESHUP plots a mesh as a set of lines.  The lines have been divided
C   --into parts as follows:
C   --   -1) Mesh boundary
C   --    0) Element block boundary
C   --    n) Interior within element block 'n'
C   --
C   --The line type and the display options are used to select the line color
C   --and the line width.
C   --
C   --For a 3D mesh with interior faces requested, the interior
C   --faces are drawn before the surface lines.
C   --
C   --Parameters:
C   --   WIDLIN - IN - true iff mesh lines should be wide versus narrow
C   --   MSHLIN - IN - the mesh lines to display (as in /MSHOPT/)
C   --   MLNTYP - IN - the line type of lines (as in /MSHOPT/)
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   LENL - IN - the cumulative line counts by element block
C   --   LINSET - IN - the sorted line set
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
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, LLNSET of /D3NUMS/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'debug.blk'

      LOGICAL WIDLIN
      INTEGER MLNTYP(-1:1)
      INTEGER IELBST(NELBLK)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      INTEGER BLKCOL(0:NELBLK)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IDELB(*)

      LOGICAL GRABRT
      LOGICAL SELONL
      LOGICAL DRAWIT
      LOGICAL LDUM
      INTEGER IDUM(1)

C   --Set mesh display parameters

      IF (MSHLIN .GE. MSHSEL) THEN
         NELB = NELBLK
      ELSE IF (MSHLIN .GE. MSHDIV) THEN
         NELB = 0
      ELSE IF (MSHLIN .GE. MSHBOR) THEN
         NELB = -1
      ELSE
         NELB = -2
      END IF
      SELONL = (MSHLIN .LE. MSHSEL)

C   --Draw the boundary determined by the hidden line removal

      IF (IS3DIM) THEN
         DO 110 IELB = NELBLK, MIN (1, NELB+1), -1

C         --If element block selected, line will be drawn later
            IF ((NELB .GE. IELB) .AND. SELONL .AND. (IELB .GE. 1)) THEN
               DRAWIT = (IELBST(IELB) .LE. 0)
            ELSE
               DRAWIT = .TRUE.
            END IF

            IF (DRAWIT) THEN

C            --Get line color and set line size
               CALL MSHCOL ('BOUNDARY', IELB, MLNTYP, .TRUE., BLKCOL,
     &            IDELB)

C            --Draw all lines marked with LINSET(3,x) < 0

               DO 100 IL = LENL(IELB-1)+1, LENL(IELB)
                  IF (LINSET(3,IL) .LT. 0) THEN
                     N1 = LINSET(1,IL)
                     IF (LINSET(3,IL) .EQ. -1) THEN
                        N2 = LINSET(2,IL)
                     ELSE
                        N2 = -LINSET(3,IL)
                     END IF
                     IF (GRABRT ()) GOTO 160
                     CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
                  END IF
  100          CONTINUE

               CALL PLTFLU
            END IF
  110    CONTINUE
      END IF

C   --Draw the selected lines

      DO 130 IELB = NELB, -1, -1

         IF (SELONL .AND. (IELB .GE. 1)) THEN
            DRAWIT = (IELBST(IELB) .GT. 0)
         ELSE
            DRAWIT = .TRUE.
         END IF

         IF (DRAWIT) THEN

C         --Get line color and set line size
            CALL MSHCOL ('LINE', IELB, MLNTYP, WIDLIN, BLKCOL, IDELB)

C         --Draw all lines of the appropriate type

            DO 120 IL = LENL(IELB-1)+1, LENL(IELB)
               IF (IS3DIM) THEN
C               --Skip invisible lines
                  IF (LINSET(3,IL) .EQ. 0) GOTO 120
               END IF
               N1 = LINSET(1,IL)
               N2 = LINSET(2,IL)
C                     write (*,*) 'line: ', n1, n2, linset(3,il)
               IF (IS3DIM) THEN
C               --Replace invisible node with partial line node
                  IF (ABS (LINSET(3,IL)) .GT. 1) N2 = ABS (LINSET(3,IL))
               END IF

               IF (GRABRT ()) GOTO 160
               CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
  120       CONTINUE

            CALL PLTFLU
         END IF
  130 CONTINUE

      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1)
     &   .and. is3dim) then
         do 150 ielb = nelblk, -1, -1
            call mshcol ('DEBUG', ielb, mlntyp, widlin, blkcol,
     &         idelb)
            do 140 il = lenl(ielb-1)+1, lenl(ielb)
               if (linset(3,il) .lt. 0) then
                  n1 = linset(1,il)
                  if (linset(3,il) .eq. -1) then
                     n2 = linset(2,il)
                  else
                     n2 = -linset(3,il)
                  end if
                  if (grabrt ()) goto 160
                  call mpd2vc (1, xn(n1), yn(n1), xn(n2), yn(n2))
               end if
  140       continue
            call pltflu
  150    continue
      end if

C   --Reset line type and size
      CALL MSHCOL ('RESET', -999, IDUM, LDUM, BLKCOL, IDELB)

      RETURN

  160 CONTINUE
C   --Reset line type and size
      CALL MSHCOL ('RESET', -999, IDUM, LDUM, BLKCOL, IDELB)
      RETURN 1
      END
