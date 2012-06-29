C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: mrkess.F,v $
C Revision 1.6  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.5  2000/07/11 16:50:23  gdsjaar
C SGI had problems with logical statment function, inlined it
C
C Revision 1.4  1997/06/26 17:51:34  caforsy
C Substituted if NeedsDoubleEscape for ifdef sun.  This define is needed on cougar, sun and ibm platforms
C
C Revision 1.3  1996/05/15 16:17:47  caforsy
C Changed Imakefile.  The cray.cf file has now removed the static heap
C and the ifdef CrayArchitecture statement is no longer needed.
C Changed the size of the arrows that mark the side sets.
C
C Revision 1.2  1994/04/08 13:25:57  gdsjaar
C Removed hash mark from comments.
C
c Revision 1.1  1994/04/07  20:05:07  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
CRevision 1.3  1993/07/08  21:57:13  gdsjaar
CReplaced SUN with sun in all ifdef statements.
C
CRevision 1.2  1990/12/14  08:53:32  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE MRKESS (LENF, NLNKF, LINKF,
     &   HIDENP, HIDEF, XN, YN, ZN, XF, YF, ZF,
     &   IELBST, IX2NP, IN2ELB, DODEAD, IDN2B, IF2EL, IE2ELB,
     &   NESSET, ISSESS, IDESS, NEESS, NNESS, IXEESS, IXNESS,
     &   LTEESS, LTNESS, *)
C=======================================================================

C   --*** MRKESS *** (SETS) Mark side set on mesh
C   --   Written by Amy Gilkey - revised 04/08/88
C   --
C   --MRKESS marks all the nodes in selected side sets on the mesh.
C   --The nodes in each side are connected in the color assigned
C   --to the side set.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the element center coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IX2NP - IN - the node number for each mesh index
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   DODEAD - IN - mark dead nodes iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   IF2EL - IN - the element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   NESSET - IN - the number of selected side sets
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTNESS - IN - the nodes for all sets
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KHCHSZ=1, KSCHSZ=2)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL HIDENP(*)
      LOGICAL HIDEF(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XF(*), YF(*), ZF(*)
      INTEGER IELBST(NELBLK)
      INTEGER IX2NP(NUMNPF)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER IF2EL(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER ISSESS(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)

      LOGICAL GRABRT
      LOGICAL PLTGTT, PLTSTT, LDUM
      LOGICAL EXISTS
      INTEGER IFACES(10)
      LOGICAL NPOK
      INTEGER NOK(10)
      LOGICAL ARROW

      EXISTS (M) = (MOD(M,2) .NE. 0)

C   --Set the symbol size, do not exit before resetting
      LDUM = PLTGTT (KSCHSZ, SZSYM)
      LDUM = PLTSTT (KSCHSZ, 1.0 * SZSYM)

      DO 150 IX = 1, NESSET
         IESS = ISSESS(IX)
         IF (NEESS(IESS) .LE. 0) THEN
            CALL PRTERR ('ERROR', 'Number of elements per side = 0')
            GOTO 150
         END IF
         NPPSID = NNESS(IESS) / NEESS(IESS)
         NPOK = (NNESS(IESS) .EQ. (NEESS(IESS) * NPPSID))
         IF (.NOT. NPOK) THEN
            CALL PRTERR ('ERROR',
     &         'Number of elements and nodes do not match')
            GOTO 150
         END IF
         IF (NPPSID .GT. 10) THEN
            CALL PRTERR ('ERROR', 'Number of nodes per side > 10')
            GOTO 150
         END IF

         CALL GRCOLR (IESS)

         IXE0 = IXEESS(IESS) - 1

         DO 140 IEE = 1, NEESS(IESS)
            IEL = LTEESS(IXE0+IEE)
            IF (IELBST(IE2ELB(IEL)) .GT. 0) THEN
               IXN0 = IXNESS(IESS) + (IEE-1) * NPPSID - 1

               DO 100 INE = 1, NPPSID
                  NOK(INE) = 0
                  INP = LTNESS(IXN0+INE)
                  IF (IS3DIM) THEN
                     IF (HIDENP(INP)) GOTO 100
                  END IF

                  if (IN2ELB(inp) .GE. 0) then
                     nok(ine) = inp
                  else if (dodead) then
                     if (idn2b(inp) .ge. 0) nok(ine) = inp
                  end if
 100           CONTINUE

               IF (NPPSID .GT. 2) THEN
                  N2 = NOK(NPPSID)
               ELSE
                  N2 = 0
               END IF

               ARROW = .TRUE.
               DO 110 INE = 1, NPPSID
                  N1 = N2
                  N2 = NOK(INE)
                  IF ((N1 .GT. 0) .AND. (N2 .GT. 0)) THEN
                     IF (GRABRT ()) GOTO 160
                     IF (ARROW) THEN
                        CALL MP2VC (1, XN(N1), YN(N1), XN(N2), YN(N2),
     &                     DX1, DY1, DX2, DY2, MASK)
                        IF (EXISTS (MASK))
     &                    CALL PLTARR (DX1, DY1, DX2, DY2, 1.0, 0.015)
                        ARROW = .FALSE.
                     ELSE
                        CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
                     END IF
C????                  ELSE IF (N2 .GT. 0) THEN
C????                     I = INE + 1
C????                     IF (I .GT. NPPSID) I = 1
C????                     IF (NOK(I) .LE. 0) THEN
C????                        CALL MPD2SY (1, XN(N2), YN(N2), '\CX')
C????                     END IF
                  END IF
  110          CONTINUE

               if (cdebug .eq. 'SSETS') then
                  CALL FNDE2F (IEL, LENF, IF2EL, NFARY, IFACES, IELB)
                  DO 130 N = 1, NFARY
                     IFAC = IFACES(N)
                     IF (IS3DIM) THEN
                        IF (HIDEF(IFAC)) GOTO 130
                     END IF
                     IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
                     DO 120 INE = 1, NPPSID
                        INP = LTNESS(IXN0+INE)
                        IF (LOCINT (INP, NLNKF(IELB), LINKF(IXL))
     &                     .LE. 0) GOTO 130
  120                CONTINUE

                     IF (GRABRT ()) GOTO 160
#if NeedsDoubleEscape
                     CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\\CX')
#else
                     CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\CX')
#endif
  130             CONTINUE
               end if
            END IF
  140    CONTINUE

         CALL PLTFLU
  150 CONTINUE

C   --Reset the symbol size, do not exit before resetting
      LDUM = PLTSTT (KSCHSZ, SZSYM)
      RETURN

  160 CONTINUE
C   --Reset the symbol size, do not exit before resetting
      LDUM = PLTSTT (KSCHSZ, SZSYM)
      RETURN 1
      END
