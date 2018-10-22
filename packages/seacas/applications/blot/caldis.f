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

C $Log: caldis.f,v $
C Revision 1.4  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.2  1996/06/21 16:06:58  caforsy
C Ran ftnchek and removed unused variables.  Reformat output for list
C var, list global, and list name.
C
C Revision 1.1  1994/04/07 19:55:09  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:47  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CALDIS (NNENUM, NENUM, XNE, YNE, ZNE, SDISTS)
C=======================================================================

C   --*** CALDIS *** (SPLOT) Calculate node/element distances
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --CALDIS calculates the distances between selected nodes or elements.
C   --
C   --Parameters:
C   --   NNENUM - IN  - the number of selected node/element numbers
C   --   NENUM  - IN  - the selected node/element numbers
C   --   XNE,
C   --   YNE,
C   --   ZNE    - IN  - the coordinates of the nodes (if nodes)
C   --                  or element centroids (if elements)
C   --   SDISTS - OUT - the calculated distances
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/

      include 'dbnums.blk'

      INTEGER NENUM(NNENUM)
      REAL XNE(*), YNE(*), ZNE(*)
      REAL SDISTS(NNENUM)

      SDISTS(1) = 0.0
      INE = NENUM(1)
      XTHIS = XNE(INE)
      YTHIS = YNE(INE)
      IF (NDIM .EQ. 3) ZTHIS = ZNE(INE)

      DO 100 NUM = 2, NNENUM
         INE = NENUM(NUM)
         XLAST = XTHIS
         XTHIS = XNE(INE)
         YLAST = YTHIS
         YTHIS = YNE(INE)
         IF (NDIM .EQ. 2) THEN
            SDISTS(NUM) = SDISTS(NUM-1) +
     &         SQRT ((XTHIS-XLAST)**2 + (YTHIS-YLAST)**2)
         ELSE IF (NDIM .EQ. 3) THEN
            ZLAST = ZTHIS
            ZTHIS = ZNE(INE)
            SDISTS(NUM) = SDISTS(NUM-1) +
     &         SQRT ((XTHIS-XLAST)**2 + (YTHIS-YLAST)**2
     &         + (ZTHIS-ZLAST)**2)
         END IF
  100 CONTINUE

      RETURN
      END
