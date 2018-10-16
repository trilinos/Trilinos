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

C $Log: setmsh.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:11:41  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.2  1990/12/14  08:57:14  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE SETMSH (IVIEW, MDEF, MNUM, MLIN, LTYP,
     &   NNPS, ISNPS, NESS, ISESS, MMOD, MTYP,
     &   ISSNPS, ISSESS)
C=======================================================================

C   --*** SETMSH *** (MESH) Set display mode and type
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --SETMSH sets the mesh display options of one or all views.  If the
C   --view is "NONE" or "EMPTY", the other options are set to default
C   --options.
C   --
C   --Parameters:
C   --   IVIEW - IN - the view to be set, 0 for all
C   --   MDEF - IN - the deformed/undeformed option to be set
C   --   MNUM - IN - the view numbering to be set
C   --   MLIN - IN - the display type for the mesh lines to be set
C   --   LTYP - IN - the mesh line type for lines
C   --   NNPS - IN - the number of selected node sets
C   --   ISNPS - IN - the indices of the selected node sets
C   --   NESS - IN - the number of selected side sets
C   --   ISESS - IN - the number of selected side sets
C   --   MMOD - IN - the display mode to be set
C   --   MTYP - IN - the display mode type to be set
C   --   ISSNPS - IN/OUT - the indices of the selected node sets
C   --   ISSESS - IN/OUT - the indices of the selected side sets
C   --
C   --Common Variables:
C   --   Uses NUMNPS, NUMESS of /DBNUMG/
C   --   Sets MSHDEF, MSHNUM, MSHLIN, MLNTYP, NNPSET, NESSET of /MSHOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'dbnumgq.blk'
      include 'mshopt.blk'

      CHARACTER*(*) MDEF, MNUM
      INTEGER LTYP(-1:1)
      INTEGER ISNPS(*)
      INTEGER ISESS(*)
      CHARACTER*(*) MMOD, MTYP
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)

      INTEGER NDEFVW, IXVW

      IF (IVIEW .GE. 1) THEN
         IF ((MDEF .EQ. 'NONE') .OR. (MDEF .EQ. 'EMPTY')) THEN
            MSHDEF(IVIEW) = MDEF
            MSHNUM(IVIEW) = 'NONE'
            MSHLIN(IVIEW) = MSHNON
            CALL INIINT (3, 0, MLNTYP(-1,IVIEW))
            NNPSET(IVIEW) = 0
            CALL INIINT (NUMNPS, 0, ISSNPS(1,IVIEW))
            NESSET(IVIEW) = 0
            CALL INIINT (NUMESS, 0, ISSESS(1,IVIEW))

            CALL SETMOD (IVIEW, 'NONE', ' ')

         ELSE
            MSHDEF(IVIEW) = MDEF
            MSHNUM(IVIEW) = MNUM
            MSHLIN(IVIEW) = MLIN
            CALL CPYINT (3, LTYP(-1), MLNTYP(-1,IVIEW))
            NNPSET(IVIEW) = NNPS
            CALL CPYINT (NNPS, ISNPS, ISSNPS(1,IVIEW))
            NESSET(IVIEW) = NESS
            CALL CPYINT (NESS, ISESS, ISSESS(1,IVIEW))

            IF (MMOD .NE. ' ') CALL SETMOD (IVIEW, MMOD, MTYP)
         END IF

      ELSE
         IF ((MDEF .EQ. 'NONE') .OR. (MDEF .EQ. 'EMPTY')) THEN
            DO 100 I = 1, 4
               MSHDEF(I) = MDEF
               MSHNUM(I) = 'NONE'
               MSHLIN(I) = MSHNON
               CALL INIINT (3, 0, MLNTYP(-1,I))
               NNPSET(I) = 0
               CALL INIINT (NUMNPS, 0, ISSNPS(1,I))
               NESSET(I) = 0
               CALL INIINT (NUMESS, 0, ISSESS(1,I))

               CALL SETMOD (I, 'NONE', ' ')
  100       CONTINUE

         ELSE
            DO 110 IVW = 1, NDEFVW (.FALSE.)
               I = IXVW (.FALSE., IVW)
               MSHDEF(I) = MDEF
               MSHNUM(I) = MNUM
               MSHLIN(I) = MLIN
               CALL CPYINT (3, LTYP(-1), MLNTYP(-1,I))
               NNPSET(I) = NNPS
               CALL CPYINT (NNPS, ISNPS, ISSNPS(1,I))
               NESSET(I) = NESS
               CALL CPYINT (NESS, ISESS, ISSESS(1,I))

               IF (MMOD .NE. ' ') CALL SETMOD (I, MMOD, MTYP)
  110       CONTINUE
         END IF
      END IF

      RETURN
      END
