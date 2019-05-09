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
      SUBROUTINE MSSURF (A, LENE, NLNKE, LINKE,
     &   LENF, NLNKF, KLINKF, KIF2EL, NAMELB)
C=======================================================================

C   --*** MSSURF *** (MESH) Convert elements into faces
C   --   Written by Amy Gilkey - revised 03/29/88
C   --
C   --MSSURF converts 3D elements into faces and packs the faces by
C   --element block.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity; connectivity all zero
C   --      if element is undefined
C   --   LENF - OUT - the cumulative face counts by element block
C   --   NLNKF - OUT - the number of nodes per face
C   --   KLINKF - OUT - the dynamic memory index of the connectivity
C   --      for all faces
C   --   KIF2EL - OUT - the dynamic memory index of the element number
C   --      of each face
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      DIMENSION A(*)
      INTEGER LENE(0:*)
      INTEGER NLNKE(*)
      INTEGER LINKE(*)
      INTEGER LENF(0:*)
      INTEGER NLNKF(*)
      CHARACTER*(MXSTLN) NAMELB(*)

      IF (.NOT. IS3DIM) THEN

         CALL CNTLNK (NELBLK, LENE, NLNKE, LENLNK, MAXFAC)
         CALL MDLONG ('LINKF', KLINKF, LENLNK)
         CALL MDLONG ('IF2EL', KIF2EL, MAXFAC)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL SURF2D (LENE, NLNKE, LINKE,
     &      LENF, NLNKF, A(KLINKF), A(KIF2EL),
     &      MAXFAC, LENLNK)

      ELSE

         CALL MDRSRV ('NLKSCR', KNLNKS, NELBLK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL CNTLK3 (NELBLK, LENE, NLNKE, LSCLNK, NSCFAC, NAMELB,
     $        A(KNLNKS))

         CALL MDRSRV ('LENSCR', KLENS, 1 + NELBLK)
         CALL MDRSRV ('LNKSCR', KLINKS, LSCLNK)
         CALL MDRSRV ('IFSCR', KF2ES, 2 * NSCFAC)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL FACE3D (A, LENE, NLNKE, LINKE, MAXFAC, LENLNK,
     $        A(KLENS), A(KNLNKS), A(KLINKS), A(KF2ES), NAMELB)

         CALL CNTLNK (NELBLK, A(KLENS), A(KNLNKS), LSCLNK, NSCFAC)
         CALL MDLONG ('LNKSCR', KLINKS, LSCLNK)
         CALL MDLONG ('IFSCR', KF2ES, 2 * NSCFAC)

         CALL MDLONG ('LINKF', KLINKF, LENLNK)
         CALL MDLONG ('IF2EL', KIF2EL, MAXFAC)
C      --MSSTEP and MSGEOM uses MDFIND to find IF2EL2
         CALL MDRSRV ('IF2EL2', KIF2E2, MAXFAC)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL SURF3D (A(KLENS), A(KNLNKS), A(KLINKS), A(KF2ES),
     &      LENF, NLNKF, A(KLINKF), A(KIF2EL), A(KIF2E2),
     &      MAXFAC, LENLNK)

         CALL MDDEL ('LENSCR')
         CALL MDDEL ('NLKSCR')
         CALL MDDEL ('LNKSCR')
         CALL MDDEL ('IFSCR')
      END IF

      CALL MDLONG ('LINKF', KLINKF, LENLNK)
      CALL MDLONG ('IF2EL', KIF2EL, MAXFAC)
      IF (IS3DIM) CALL MDLONG ('IF2EL2', KIF2E2, MAXFAC)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

  100 CONTINUE
      RETURN
      END
