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

C $Log: lnread.f,v $
C Revision 1.4  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  1998/06/12 15:53:24  gdsjaar
C 1. Problem with TIMES array. Blot accesses a dummy timestep even if
C there were no timesteps on the database. Array wasn't allocated, so
C writing off into never-never land.
C
C 2. Inconsistency among some common blocks. Some places weren't using
C the include but had the definition hardwired in. Removed those.
C
C 3. Added 'EXTERNAL BLKDAT' to all routines that used data values set
C in BLKDAT
C
C 4. Cleanup of some A vs. IA argument passing.
C
C Revision 1.2  1994/04/08 13:25:56  gdsjaar
C Removed hash mark from comments.
C
c Revision 1.1  1994/04/07  20:04:30  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
c Revision 1.2  1990/12/14  08:53:06  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE LNREAD (A, NPTIMS, NPTIMW, IPTIMS, TIMES, WHOTIM,
     &   XLN, YLN, ZLN)
C=======================================================================

C   --*** LNREAD *** (PATHLN) Read pathline data from database
C   --   Written by Amy Gilkey - revised 05/27/88
C   --
C   --LNREAD reads the database and stores the pathline data.
C   --
C   --This routine manipulates dynamic memory, so check after return.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NPTIMS - IN - the number of selected steps (for history variables)
C   --   NPTIMW - IN - the number of whole selected steps
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the time step times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   XLN, YLN, ZLN - OUT - the pathline data
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses NLNCRV, ILVID of /LNVARS/

      include 'dbnums.blk'
      include 'lnvars.blk'

      DIMENSION A(*)
      INTEGER IPTIMS(NPTIMS)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      REAL XLN(NPTIMS,NLNCRV), YLN(NPTIMS,NLNCRV), ZLN(NPTIMS,NLNCRV)

      LOGICAL NEEDHV, NEEDGV, NEEDNV, NEEDEV
      CHARACTER TYP

C   --Determine which types of variables are needed

      NEEDHV = .FALSE.
      NEEDGV = .FALSE.
      NEEDNV = .FALSE.
      NEEDEV = .FALSE.
      LDATA = 0

      DO 100 NP = 1, NLNCRV
         CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
         IF (TYP .EQ. 'H') THEN
            NEEDHV = .TRUE.
            LDATA = MAX (LDATA, NVARHI)
         ELSE IF (TYP .EQ. 'G') THEN
            NEEDGV = .TRUE.
            LDATA = MAX (LDATA, NVARGL)
         ELSE IF (TYP .EQ. 'N') THEN
            NEEDNV = .TRUE.
            LDATA = MAX (LDATA, NUMNP)
         ELSE IF (TYP .EQ. 'E') THEN
            NEEDEV = .TRUE.
            LDATA = MAX (LDATA, NUMEL)
         END IF
  100 CONTINUE

C   --Reserve memory for data record

      CALL MDRSRV ('DATA', KDATA, LDATA)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

C   --Transfer element variables onto random file (for efficiency)

      MXSTEP = 0
      DO 110 NPT = 1, NPTIMS
         ISTEP = IPTIMS(NPT)
         IF (WHOTIM(ISTEP)) MXSTEP = ISTEP
  110 CONTINUE

      IF (NEEDEV .AND. (MXSTEP .GT. 0)) THEN
C????         CALL LNTRND (A, MXSTEP, 'E', NUMEL, NVAREL, A(KDATA))
      END IF

      NPTIMW = 0
      DO 120 NPT = 1, NPTIMS

         ISTEP = IPTIMS(NPT)
         IF (WHOTIM(NPT)) NPTIMW = NPTIMW + 1

C      --Read and store variable data to be plotted

         IF (NEEDHV) THEN
            CALL LNSTOR (A, ISTEP, 'H', NVARHI, 1,
     &         NPT, NPTIMS, XLN, YLN, ZLN, A(KDATA))
         END IF

         IF (WHOTIM(ISTEP)) THEN
            IF (NEEDGV) THEN
               CALL LNSTOR (A, ISTEP, 'G', NVARGL, 1,
     &            NPTIMW, NPTIMS, XLN, YLN, ZLN, A(KDATA))
            END IF

            IF (NEEDNV) THEN
               CALL LNSTOR (A, ISTEP, 'N', NUMNP, NVARNP,
     &            NPTIMW, NPTIMS, XLN, YLN, ZLN, A(KDATA))
            END IF

            IF (NEEDEV) THEN
               CALL LNSTOR (A, ISTEP, 'E', NUMEL, NVAREL,
     &            NPTIMW, NPTIMS, XLN, YLN, ZLN, A(KDATA))
            END IF
         END IF

  120 CONTINUE

      CALL MDDEL ('DATA')

  130 CONTINUE
      RETURN
      END
