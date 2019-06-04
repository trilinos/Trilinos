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

C $Log: cpymsh.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:57:21  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.2  1990/12/14  08:49:06  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE CPYMSH (INEW, IOLD, ISSNPS, ISSESS)
C=======================================================================

C   --*** CPYMSH *** (MESH) Copy mesh display parameters from another view
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --CPYMSH copies the mesh display options from one view to another
C   --view or sets the view to undefined.
C   --
C   --Parameters:
C   --   INEW - IN - the view to be set
C   --   IOLD - IN - the view to be copied; <=0 to set to undefined
C   --   ISSNPS - IN/OUT - the indices of the selected node sets
C   --   ISSESS - IN/OUT - the indices of the selected side sets
C   --
C   --Common Variables:
C   --   Uses NUMNPS, NUMESS of /DBNUMG/
C   --   Sets MSHDEF, MSHNUM, MSHLIN, MLNTYP, NNPSET, NESSET of /MSHOPT/
C   --   Uses MODDET, MODTYP of /DETOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'dbnumgq.blk'
      include 'mshopt.blk'
      include 'detopt.blk'

      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)

      CHARACTER CDUM

      IF (IOLD .GE. 1) THEN
         CALL SETMSH (INEW,
     &      MSHDEF(IOLD), MSHNUM(IOLD), MSHLIN(IOLD), MLNTYP(-1,IOLD),
     &      NNPSET(IOLD), ISSNPS(1,IOLD), NESSET(IOLD), ISSESS(1,IOLD),
     &      MODDET(IOLD), MODTYP(IOLD), ISSNPS, ISSESS)

      ELSE
         CALL SETMSH (INEW,
     &      'NONE', CDUM, IDUM, IDUM,
     &      IDUM, IDUM, IDUM, IDUM,
     &      'NONE', CDUM, ISSNPS, ISSESS)
      END IF

      RETURN
      END
