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

C $Log: allal3.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 19:54:37  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:32  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ALLAL3 (LENF, IF2EL, IF2EL2, IE2ELB, NEWELB)
C=======================================================================

C   --*** ALLAL3 *** (MESH) Adjust for elements all alive (3D)
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --ALLAL3 adjusts the face array to reflect all elements alive.
C   --
C   --An element birth causes all the faces of the element to be moved
C   --as follows:
C   --   Dead    -> surface (if one alive element for face) or
C   --              interior (if two alive elements)
C   --   Surface -> interior
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   IF2EL2 - IN - the secondary element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   NEWELB - OUT - size = LENF(NELBLK+1)
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK+2)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER NEWELB(*)

C   --Change half-dead surface faces back to interior faces

      DO 110 IELB = 1, NELBLK
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IF2EL2(IFAC) .GT. 0) THEN
               NEWELB(IFAC) = NELBLK+1
            END IF
  100    CONTINUE
  110 CONTINUE

C   --Leave interior faces alone

      IELB = NELBLK+1
      DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
         NEWELB(IFAC) = IELB
  120 CONTINUE

C   --Change dead faces back to surface or interior faces

      IELB = NELBLK+2
      DO 130 IFAC = LENF(IELB-1)+1, LENF(IELB)
         IF (IF2EL2(IFAC) .LE. 0) THEN
            NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
         ELSE
            NEWELB(IFAC) = NELBLK+1
         END IF
  130 CONTINUE

      RETURN
      END
