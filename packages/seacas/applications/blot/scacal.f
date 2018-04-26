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

C $Log: scacal.f,v $
C Revision 1.3  2009/03/25 12:36:47  gdsjaar
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
C Revision 1.1  1994/04/07 20:10:30  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:46  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCACAL (NAME, IVAR, USESEL, IELBST,
     &   ISTMN, ICALC)
C=======================================================================

C   --*** SCACAL *** (BLOT) Return variable scaling if already calculated
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --SCACAL determines if the variables has already been scaled.  If not,
C   --a scaling message is printed.  For element variables, the message
C   --is printed if all element blocks are not selected.
C   --
C   --Parameters:
C   --   NAME - IN - the variable name
C   --   IVAR - IN - the variable index
C   --   USESEL - IN - use the element blocks selected array iff true,
C   --      else all selected
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISTMN - IN - the "state" of the minimums and maximums:
C   --      <0 if not calculated
C   --       0 if no elements in element block
C   --      +n if calculated
C   --   ICALC - OUT - the scaling flag
C   --      0 = all element blocks selected (if element) and min/max calculated
C   --      1 = not all element blocks selected (if element), but min/max
C   --          calculated for each selected element block
C   --      2 = min/max not calculated for some element block
C   --
C   --Common Variables:
C   --   Uses NELBLK, NVARNP, NVAREL of /DBNUMS/

      include 'dbnums.blk'

      CHARACTER*(*) NAME
      LOGICAL USESEL
      INTEGER IELBST(*)
      INTEGER ISTMN(0:*)

      CHARACTER TYP

C   --Get variable type

      CALL DBVTYP_BL (IVAR, TYP, IDUM)

      IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')
     &   .OR. (TYP .EQ. 'N')) THEN

C      --History, global or nodal variable, check min/max already computed

         IF (ISTMN(0) .GE. 0) THEN
            ICALC = 0
         ELSE
            ICALC = 2
         END IF

      ELSE IF (TYP .EQ. 'E') THEN

C      --Element variable, check if all element blocks selected and if
C      --all-blocks min/max already computed

         IF (USESEL) THEN
            CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
         ELSE
            NUMSEL = NELBLK
         END IF

         IF ((NUMSEL .GE. NELBLK) .AND. (ISTMN(0) .GE. 0)) THEN
            ICALC = 0

         ELSE

C         --If not all element blocks selected, determine if each min/max for
C         --variable already computed

            ICALC = 1
            DO 100 IELB = 1, NELBLK
               IF (ISTMN(IELB) .LT. 0) ICALC = 2
  100       CONTINUE
         END IF

      END IF

      IF (ICALC .GT. 0) THEN
         WRITE (*, 10000) 'Scaling variable ', NAME(:LENSTR(NAME))
      END IF

      RETURN
10000  FORMAT (1X, 5A)
      END
