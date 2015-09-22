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

C $Log: sorlnk.f,v $
C Revision 1.3  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:13:50  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:05  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SORLNK (A, NSETS, NEWELB,
     &   LENF, NLNKF, LINKF, IF2EL, IF2EL2, IE2ELB)
C=======================================================================

C   --*** SORLNK *** (MESH) Sort faces
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --SORLNK resorts all the faces into new element block sets.  This
C   --is done by coping the faces, then renaming the arrays.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   LINKF - for length only
C   --   IF2EL - for length only
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NSETS - IN - the number of element block sets in LINKF
C   --   NEWELB - IN - the new element block set for each face
C   --   LENF - IN/OUT - the cumulative face counts by element block for LINKF
C   --   NLNKF - IN - the number of nodes for each element block
C   --   LINKF - IN/OUT - the unsorted connectivity for all faces
C   --   IF2EL - IN/OUT - the element number of each face in LINKF
C   --   IF2EL2 - IN/OUT - the secondary element number of each face in LINKF
C   --   IE2ELB - IN - the element for each element block
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      DIMENSION A(*)
      INTEGER NEWELB(*)
      INTEGER LENF(0:NSETS)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(*)

      LOGICAL SOMMOV

C   --Check that some faces will move

      SOMMOV = .FALSE.
      DO 110 IELB = 1, NSETS
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IELB .NE. IABS (NEWELB(IFAC))) THEN
               SOMMOV = .TRUE.
               GOTO 120
            END IF
  100    CONTINUE
  110 CONTINUE
  120 CONTINUE

      IF (SOMMOV) THEN
         CALL MDRSRV ('XLENF', KLENF, 1+NSETS)
         CALL MDFIND ('LINKF', IDUM, LENLNK)
         CALL MDRSRV ('XLINKF', KLINKF, LENLNK)
         CALL MDFIND ('IF2EL', IDUM, LENF2E)
         CALL MDRSRV ('XIF2EL', KIF2EL, LENF2E)
         IF (IS3DIM) CALL MDRSRV ('XIF2E2', KIF2E2, LENF2E)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 130

         CALL CPYINT (1+NSETS, LENF, A(KLENF))
         CALL CPYINT (LENLNK, LINKF, A(KLINKF))
         CALL CPYINT (LENF2E, IF2EL, A(KIF2EL))
         IF (IS3DIM) CALL CPYINT (LENF2E, IF2EL2, A(KIF2E2))

         CALL SORLNX (NSETS, NEWELB, IE2ELB,
     &      NLNKF, A(KLENF), A(KLINKF), A(KIF2EL), A(KIF2E2),
     &      LENF, LINKF, IF2EL, IF2EL2)

         CALL MDDEL ('XLENF')
         CALL MDDEL ('XLINKF')
         CALL MDDEL ('XIF2EL')
         IF (IS3DIM) CALL MDDEL ('XIF2E2')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 130
      END IF
  130 CONTINUE
      RETURN
      END
