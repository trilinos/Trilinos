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

C $Log: cntvar.f,v $
C Revision 1.3  2009/03/25 12:36:43  gdsjaar
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
C Revision 1.1  1994/04/07 19:56:57  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:52  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CNTVAR (MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR)
C=======================================================================

C   --*** CNTVAR *** (DETOUR) Count the number of variables needed
C   --   Written by Amy Gilkey - revised 03/03/88
C   --
C   --CNTVAR counts the number of nodal and element database variables
C   --needed for the display modes.
C   --
C   --Parameters:
C   --   MODDET - IN - the modes for all views (as in /DETOPT/)
C   --   MODTYP - IN - the mode types for all views (as in /DETOPT/)
C   --   IDTVAR - IN - the current variables
C   --   NNDVAR - OUT - the number of nodal variables needed
C   --   NEDVAR - OUT - the number of element variables needed
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/

      include 'dbnums.blk'

      CHARACTER*(*) MODDET(4), MODTYP(4)
      INTEGER IDTVAR(4)

      INTEGER NDEFVW, IXVW
      CHARACTER TYP

      NNDVAR = 0
      NEDVAR = 0
      DO 100 IVW = 1, NDEFVW (.FALSE.)
         IVIEW = IXVW (.FALSE., IVW)
         IF (MODDET(IVIEW) .EQ. 'CONTOUR') THEN
            NNDVAR = MAX (NNDVAR, 1)
            CALL DBVTYP_BL (IDTVAR(1), TYP, IDUM)
            IF (TYP .EQ. 'E') NEDVAR = MAX (NEDVAR, 1)
         ELSE IF (MODDET(IVIEW) .EQ. 'ELEMCONT') THEN
            NEDVAR = MAX (NEDVAR, 1)
         ELSE IF (MODDET(IVIEW) .EQ. 'VECTOR') THEN
            IF (MODTYP(IVIEW) .EQ. 'NODE') THEN
               NNDVAR = MAX (NNDVAR, NDIM)
            ELSE IF (MODTYP(IVIEW) .EQ. 'ELEMENT') THEN
               NEDVAR = MAX (NEDVAR, NDIM)
            ELSE IF ((MODTYP(IVIEW) .EQ. 'SIGMAX')
     &         .OR. (MODTYP(IVIEW) .EQ. 'SIGMIN')) THEN
               NEDVAR = MAX (NEDVAR, 3)
            END IF
         ELSE IF (MODDET(IVIEW) .EQ. 'SYMBOL') THEN
            NEDVAR = MAX (NEDVAR, 1)
         ELSE IF (MODDET(IVIEW) .EQ. 'GAUSS') THEN
            NEDVAR = MAX (NEDVAR, 4)
         END IF
  100 CONTINUE

      RETURN
      END
