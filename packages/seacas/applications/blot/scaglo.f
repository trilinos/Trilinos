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

C $Log: scaglo.f,v $
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
C Revision 1.1  1994/04/07 20:10:39  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:52  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCAGLO (A, VAR, WHOTIM,
     &   VALMIN, ISTMIN, VALMAX, ISTMAX)
C=======================================================================

C   --*** SCAGLO *** (BLOT) Scale all global variables
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAGLO reads the values for the global variables from the database
C   --and finds the minimum and maximum values.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index (for GETVAR)
C   --   VAR - SCRATCH - the variable array
C   --   WHOTIM - IN - true iff time step is a whole (versus history) time step
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value for each variable
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum
C   --      value for each variable
C   --
C   --Common Variables:
C   --   Uses NVARGL, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      REAL VAR(NVARGL)
      LOGICAL WHOTIM(*)
      REAL VALMIN(NVARGL), VALMAX(NVARGL)
      INTEGER ISTMIN(NVARGL), ISTMAX(NVARGL)

      CALL DBVIX_BL ('G', 1, IVAR)

      DO 110 ISTEP = 1, NSTEPS
         IF (.NOT. WHOTIM(ISTEP)) GOTO 110

C      --Read the variables

         CALL GETVAR (A, IVAR, -999, ISTEP, NVARGL, VAR)

C      --Find minimum and maximum variable values for variable

         DO 100 IXVAR = 1, NVARGL

            IF (ISTEP .EQ. 1) THEN
               VALMIN(IXVAR) = VAR(IXVAR)
               ISTMIN(IXVAR) = ISTEP
               VALMAX(IXVAR) = VAR(IXVAR)
               ISTMAX(IXVAR) = ISTEP
            ELSE IF (VALMIN(IXVAR) .GT. VAR(IXVAR)) THEN
               VALMIN(IXVAR) = VAR(IXVAR)
               ISTMIN(IXVAR) = ISTEP
            ELSE IF (VALMAX(IXVAR) .LT. VAR(IXVAR)) THEN
               VALMAX(IXVAR) = VAR(IXVAR)
               ISTMAX(IXVAR) = ISTEP
            END IF

  100    CONTINUE
  110 CONTINUE

      RETURN
      END
