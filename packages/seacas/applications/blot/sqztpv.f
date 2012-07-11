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

C $Log: sqztpv.f,v $
C Revision 1.3  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1999/03/09 19:41:43  gdsjaar
C Fixed missing parameter definition of MSHBOR in blotII2.f
C
C Cleaned up parameters and common blocks in other routines.
C
C Revision 1.1  1994/04/07 20:15:43  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.2  1990/12/14  08:58:44  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE SQZTPV (NPTIMS, IPTIMS, WHOTIM, NPTS, PLTVAL)
C=======================================================================

C   --*** SQZTPV *** (TPLOT) Compress curves for whole time steps
C   --   Written by Amy Gilkey - revised 11/11/87
C   --
C   --SQZTPV compresses the values for TPLOT curves so that only the values
C   --for whole time steps are in the curve.  Curves of history variables
C   --only are left intact.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   NPTS - OUT - the number of points on each curve
C   --   PLTVAL - IN/OUT - the plot data;
C   --      PLTVAL(x,NTPVAR+1) holds the times if TIMPLT
C   --      PLTVAL(x,NTPVAR+2) holds the compressed times if TIMPLT and needed
C   --
C   --Common Variables:
C   --   Uses NTPCRV, NTPVAR, TIMPLT of /TPVARS/

      include 'params.blk'
      include 'tpvars.blk'

      INTEGER IPTIMS(*)
      LOGICAL WHOTIM(*)
      INTEGER NPTS(NTPVAR+2)
      REAL PLTVAL(NPTIMS,NTPVAR+2)

      CHARACTER TYPX, TYPY

      NSQZ = 0
      NPTIMW = NWHSEL (NPTIMS, IPTIMS, WHOTIM)

      IF (TIMPLT) TYPX = 'H'
      N = 0
      DO 100 NP = 1, NTPCRV
         IF (.NOT. TIMPLT) THEN
            N = N + 1
            CALL DBVTYP_BL (ITVID(N), TYPX, IDUM)
         END IF
         N = N + 1
         CALL DBVTYP_BL (ITVID(N), TYPY, IDUM)
         IF ((NPTIMS .NE. NPTIMW) .AND.
     &      ((TYPX .NE. 'H') .OR. (TYPY .NE. 'H'))) THEN
            NSQZ = NSQZ + 1
            IF (.NOT. TIMPLT) THEN
               CALL SQZLGV (NPTIMS, IPTIMS, WHOTIM,
     &            PLTVAL(1,N-1), NPTS(N-1), PLTVAL(1,N-1))
            END IF
            CALL SQZLGV (NPTIMS, IPTIMS, WHOTIM,
     &         PLTVAL(1,N), NPTS(N), PLTVAL(1,N))
         ELSE
            IF (.NOT. TIMPLT) THEN
               NPTS(N-1) = NPTIMS
            END IF
            NPTS(N) = NPTIMS
         END IF
  100 CONTINUE

      IF ((NSQZ .GT. 0) .AND. TIMPLT) THEN
         CALL SQZLGV (NPTIMS, IPTIMS, WHOTIM,
     &      PLTVAL(1,NTPVAR+1), IDUM, PLTVAL(1,NTPVAR+2))
      END IF

      RETURN
      END
