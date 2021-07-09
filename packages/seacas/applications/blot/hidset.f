C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDSET (LINSET, LENL, HIDENP, ALLHOL, IELBST,
     &   IEDSET, NEDGES, IPSET, NPART)
C=======================================================================

C   --*** HIDSET *** (MESH) Identify partial line set
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --HIDSET identifies those lines which connect visible to hidden nodes
C   --overlapped by a face (versus hidden nodes with the normal
C   --pointing away).  The hidden nodes are moved to LINSET(2,x) and a
C   --list of indices into the LINSET array is made up.
C   --
C   --Parameters:
C   --   LINSET - IN/OUT - the sorted line set
C   --   LENL - IN - the length of the line set, for each element block
C   --   HIDENP - IN - node status (as in HIDDEN)
C   --   ALLHOL - IN - true iff mesh display is hollow
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IEDSET - IN - the edge line set; (0) = face defining edge
C   --   NEDGES - IN - the number of edges
C   --   IPSET - OUT - the indices of the partial line set
C   --   NPART - OUT - the number of lines in the partial line set
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses LLNSET of /D3NUMS/

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENL(-2:NELBLK)
      INTEGER LINSET(LLNSET,*)
      INTEGER HIDENP(*)
      LOGICAL ALLHOL
      INTEGER IELBST(NELBLK)
      INTEGER IEDSET(0:2,*)
      INTEGER IPSET(*)

      LOGICAL EDGONL

C   --Get all line sets that are partially hidden, put the hidden node
C   --in LINSET(2,x)

      NPART = 0
      DO 140 IELB = -1, NELBLK

C      --Flag whether partial lines of this element block must be in the
C      --edge set

         EDGONL = .FALSE.
         IF (IELB .GE. 1) THEN
            IF (ALLHOL .OR. (IELBST(IELB) .LE. 0)) EDGONL = .TRUE.
         END IF

         DO 130 IL = LENL(IELB-1)+1, LENL(IELB)
            IF ((HIDENP(LINSET(1,IL)) .LE. KNVIS) .OR.
     &         (HIDENP(LINSET(2,IL)) .LE. KNVIS)) THEN
               IF ((HIDENP(LINSET(1,IL)) .EQ. KNFOVR) .OR.
     &            (HIDENP(LINSET(2,IL)) .EQ. KNFOVR)) THEN

C               --Get hidden node in LINSET(2,x)
                  IF (HIDENP(LINSET(1,IL)) .GT. KNVIS) THEN
                     I = LINSET(1,IL)
                     LINSET(1,IL) = LINSET(2,IL)
                     LINSET(2,IL) = I
                  END IF

                  IF (EDGONL) THEN
C                  --Eliminate the partial line if not in the edge set
                     IV = LINSET(1,IL)
                     IH = LINSET(2,IL)
                     DO 100 IEDG = 1, NEDGES
                        IF (IV .EQ. IEDSET(1,IEDG)) THEN
                           IF (IH .EQ. IEDSET(2,IEDG)) GOTO 110
                        END IF
  100                CONTINUE
                     GOTO 120
  110                CONTINUE
                  END IF

C               --Add the line to the partial line set
                  NPART = NPART + 1
                  IPSET(NPART) = IL
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
  140 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'partial line set =', npart

      RETURN
      END
