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

C $Log: fnde2f.f,v $
C Revision 1.3  2009/03/25 12:36:44  gdsjaar
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
C Revision 1.1  1994/04/07 20:01:08  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:45  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE FNDE2F (IEL, LENF, IF2EL, NQARY, IFACES, IELB)
C=======================================================================

C   --*** FNDE2F *** (MESH) Find faces that make up element
C   --   Written by Amy Gilkey - revised 07/06/87
C   --
C   --FNDE2F finds the surface faces that make up the given element.
C   --
C   --Parameters:
C   --   IEL - IN - the element number
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   NQARY - IN/OUT - input as the maximum length of the IFACES array;
C   --      output as the length of the IFACES array
C   --   IFACES - OUT - the surface faces that make up the element
C   --   IELB - OUT - the element block of the element quarilaterals
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/

      include 'dbnums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER IF2EL(*)
      INTEGER IFACES(*)

      NQARY = 0
      DO 110 IELB = 1, NELBLK
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IEL .EQ. IF2EL(IFAC)) THEN
               NQARY = NQARY + 1
               IFACES(NQARY) = IFAC
            END IF
  100    CONTINUE
         IF (NQARY .GT. 0) GOTO 120
  110 CONTINUE
      IELB = 1

  120 CONTINUE
      RETURN
      END
