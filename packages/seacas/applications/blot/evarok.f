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

C $Log: evarok.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:00:29  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:53  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE EVAROK (NVARS, NVAR, NELBLK, IELBST, ISEVOK, ISVOK)
C=======================================================================

C   --*** EVAROK *** (MESH) Get the multi-variable truth table
C   --   Written by Amy Gilkey - revised 10/29/87
C   --
C   --EVAROK creates the multi-variable truth table.  It uses the element
C   --variable truth table and the selected variables to create a table
C   --for only the selected variables.  If no element variables are given,
C   --the truth table is all true.
C   --
C   --Parameters:
C   --   NVARS - IN - the number of variables
C   --   NVAR - IN - the variable numbers, if any
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   ISVOK - OUT - the variable truth table; true iff all variables
C   --      of block j exist and are selected

      INTEGER NVAR(*)
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      LOGICAL ISVOK(NELBLK)

      CHARACTER TYP

      DO 100 IELB = 1, NELBLK
         ISVOK(IELB) = (IELBST(IELB) .GT. 0)
  100 CONTINUE

      DO 120 IVAR = 1, NVARS
         CALL DBVTYP_BL (NVAR(IVAR), TYP, ID)
         IF (TYP .EQ. 'E') THEN
            DO 110 IELB = 1, NELBLK
               IF (.NOT. ISEVOK(IELB,ID)) THEN
                  ISVOK(IELB) = .FALSE.
               END IF
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
