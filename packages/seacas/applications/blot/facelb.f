C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C=======================================================================
      SUBROUTINE FACELB (IELB, LENE, NLNKE, LINKE,
     &  NLNKSC, LINKSC, IF2ESC, MAXNPF, NPFS, NFACES, NAMELB)
C=======================================================================

C   --*** FACELB *** (MESH) Match element block faces
C   --   Written by Amy Gilkey - revised 02/24/88
C   --              Sam Key, 06/01/85
C   --
C   --FACELB makes up a list of faces for all faces in an element block.
C   --All matching faces within the block are combined.
C   --
C   --Parameters:
C   --   IELB - IN - the element block number
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the connectivity for all elements
C   --   NLNKSC - IN - the number of nodes per face
C   --   LINKSC - IN/OUT - the connectivity for all faces (this block only)
C   --   IF2ESC - IN/OUT - the element number(s) of each face in LINKSC
C   --   MAXNPF - IN - the maximum length of the NPFS entry
C   --   NPFS - SCRATCH - the list of unmatched faces containing a node;
C   --      (0,i) = the length of the list
C   --   NFACES - OUT - the number of unique faces in this element block
C   --
C   --Common Variables:
C   --   Uses NUMNP of /DBNUMS/



      include 'params.blk'
      include 'dbnums.blk'
      include 'minmax.blk'

      INTEGER LENE(0:NELBLK), LINKE(NLNKE,*)
      INTEGER LINKSC(*)
      INTEGER IF2ESC(2,*)
      INTEGER NPFS(NUMNP,0:MAXNPF)

      INTEGER LINKF1(4)
      CHARACTER*(MXSTLN) NAMELB

      LOGICAL ISTET

      istet = (namelb(:3) .eq. 'TET')

      CALL ININPF (minnod, maxnod, MAXNPF, NPFS)
      maxnod = 1
      minnod = numnp

      NOVER = 0
      NFACES = 0
      IXF = 1
      IF (NLNKE .EQ. 8) THEN
        DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
          IXEL = IEL - LENE(IELB-1)
          IF (LINKE(1,IXEL) .NE. 0) THEN

            DO 100 IFACE = 1, 6
              CALL FNODES (IFACE, LINKE(1,IXEL), LINKF1)
              IF (    (NPFS(LINKF1(1),0) .GT. 0)
     &          .AND. (NPFS(LINKF1(2),0) .GT. 0)
     &          .AND. (NPFS(LINKF1(3),0) .GT. 0)
     &          .AND. (NPFS(LINKF1(4),0) .GT. 0)) THEN
                IMATCH = MATFAC (LINKF1, MAXNPF, NPFS, iel, numnp, IERR)
              else
                imatch = 0
              endif

              IF (IMATCH .LE. 0) THEN
                NFACES = NFACES + 1
                CALL CPYINT (NLNKSC, LINKF1, LINKSC(IXF))
                IXF = IXF + NLNKSC
                IF2ESC(1,NFACES) = IEL
                IF2ESC(2,NFACES) = 0

                CALL FILNPF (NLNKSC, LINKF1, NFACES,
     &            MAXNPF, NPFS, NOVER, NUMNP)

              ELSE
                IF2ESC(2,IMATCH) = IEL
              END IF
 100        CONTINUE
          end if
 110    continue

      ELSE IF (NLNKE .eq. 4 .and. istet) THEN
        DO 210 IEL = LENE(IELB-1)+1, LENE(IELB)
          IXEL = IEL - LENE(IELB-1)

          IF (LINKE(1,IXEL) .NE. 0) THEN

            DO 200 IFACE = 1, 4
              CALL TNODES (IFACE, LINKE(1,IXEL), LINKF1)
              IF (    (NPFS(LINKF1(1),0) .GT. 0)
     &          .AND. (NPFS(LINKF1(2),0) .GT. 0)
     &          .AND. (NPFS(LINKF1(3),0) .GT. 0)) THEN
                IMATCH = MATFAT (LINKF1, MAXNPF, NPFS, iel, numnp, IERR)
              else
                imatch = 0
              endif

              IF (IMATCH .LE. 0) THEN
                NFACES = NFACES + 1
c                     CALL CPYINT (4, LINKF1, LINKSC(IXF))
                linksc(ixf+0) = linkf1(1)
                linksc(ixf+1) = linkf1(2)
                linksc(ixf+2) = linkf1(3)
                linksc(ixf+3) = linkf1(4)
                IXF = IXF + 4
                IF2ESC(1,NFACES) = IEL
                IF2ESC(2,NFACES) = 0

                CALL FILNPF (4, LINKF1, NFACES,
     &            MAXNPF, NPFS, NOVER, NUMNP)

              ELSE
                IF2ESC(2,IMATCH) = IEL
              END IF
 200        CONTINUE
          end if
 210    continue

      ELSE IF (NLNKE .LT. 8) THEN
        DO 310 IEL = LENE(IELB-1)+1, LENE(IELB)
          IXEL = IEL - LENE(IELB-1)

          IF (LINKE(1,IXEL) .NE. 0) THEN

            NFACES = NFACES+1
            CALL CPYINT (NLNKSC, LINKE(1,IXEL), LINKSC(IXF))
            IXF = IXF + NLNKSC
            IF2ESC(1,NFACES) = IEL
            IF2ESC(2,NFACES) = 0
          end if
 310    continue

      END IF

      IF (NOVER .GT. 0) THEN
        WRITE (*, 10000) 'in FACELB, MAXNPF =', MAXNPF, ', # =', NOVER
10000   FORMAT (1X, 'PROGRAM ERROR - ', A, I5, A, I5)
      END IF

      RETURN
      END
