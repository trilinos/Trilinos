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

C $Log: scanod.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:11:07  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:57:03  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCANOD (A, IVAR, VAR, WHOTIM, XN, YN, ZN,
     &   VALMIN, NUMMIN, XYZMIN, ISTMIN, VALMAX, NUMMAX, XYZMAX, ISTMAX)
C=======================================================================

C   --*** SCANOD *** (BLOT) Scale nodal variable
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCANOD reads the values for the nodal variable from the database
C   --and finds the minimum and maximum value.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index (for GETVAR)
C   --   VAR - SCRATCH - the variable array
C   --   WHOTIM - IN - true iff time step is a whole (versus history) time step
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value
C   --   NUMMIN, NUMMAX - OUT - the node number of the minimum and maximum value
C   --   XYZMIN, XYZMAX - OUT - the coordinates of NUMMIN, NUMMAX
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum value
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NVARNP, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      REAL VAR(*)
      LOGICAL WHOTIM(*)
      REAL XN(*), YN(*), ZN(*)
      REAL VALMIN, VALMAX
      INTEGER NUMMIN, NUMMAX
      REAL XYZMIN(3), XYZMAX(3)
      INTEGER ISTMIN, ISTMAX

      DO 110 ISTEP = 1, NSTEPS
         IF (.NOT. WHOTIM(ISTEP)) GOTO 110

C      --Read the variables

         CALL GETVAR (A, IVAR, -999, ISTEP, NUMNP, VAR)

C      --Find minimum and maximum variable values for variable

         IF (ISTEP .EQ. 1) THEN
            INP = 1
            VALMIN = VAR(INP)
            NUMMIN = INP
            ISTMIN = ISTEP
            VALMAX = VAR(INP)
            NUMMAX = INP
            ISTMAX = ISTEP
         END IF

         DO 100 INP = 1, NUMNP
            IF (VALMIN .GT. VAR(INP)) THEN
               VALMIN = VAR(INP)
               NUMMIN = INP
               ISTMIN = ISTEP
            ELSE IF (VALMAX .LT. VAR(INP)) THEN
               VALMAX = VAR(INP)
               NUMMAX = INP
               ISTMAX = ISTEP
            END IF
  100    CONTINUE
  110 CONTINUE

      DO 120 I = 1, 3
         XYZMIN(I) = 0.0
  120 CONTINUE
      IF (NDIM .GE. 1) XYZMIN(1) = XN(NUMMIN)
      IF (NDIM .GE. 2) XYZMIN(2) = YN(NUMMIN)
      IF (NDIM .GE. 3) XYZMIN(3) = ZN(NUMMIN)
      DO 130 I = 1, 3
         XYZMAX(I) = 0.0
  130 CONTINUE
      IF (NDIM .GE. 1) XYZMAX(1) = XN(NUMMAX)
      IF (NDIM .GE. 2) XYZMAX(2) = YN(NUMMAX)
      IF (NDIM .GE. 3) XYZMAX(3) = ZN(NUMMAX)

      RETURN
      END
