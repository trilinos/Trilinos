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

C $Log: lnstor.f,v $
C Revision 1.3  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1998/06/12 15:53:27  gdsjaar
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
C Revision 1.1  1994/04/07 20:04:35  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:53:09  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE LNSTOR (A, ISTEP, TYP, NWRDS, NVAR, NPT, NPTS,
     &   XLN, YLN, ZLN, DATA)
C=======================================================================

C   --*** LNSTOR *** (PATHLN) Read and store pathline data from database
C   --   Written by Amy Gilkey - revised 05/27/88
C   --
C   --LNSTOR reads variables from the database and stores any that are
C   --pathline data in the appropriate location.  It reads only one
C   --variable type (history, global, nodal, element) per call.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the time step number
C   --   TYP - IN - the type of variable:
C   --      'H'istory, 'G'lobal, 'N'odal, 'E'lement
C   --   NWRDS - IN - the number of words in a data record
C   --   NVAR - IN - the number of records to read
C   --   NPT - IN - the XLN, YLN, ZLN time index to fill
C   --   NPTS - IN - the maximum XLN, YLN, ZLN time index
C   --   XLN, YLN, ZLN - IN/OUT - the pathline data array
C   --   DATA - SCRATCH - size = NWRDS
C   --
C   --Common Variables:
C   --   Uses NLNCRV, ILVNE, ILVID of /LNVARS/

      include 'lnvars.blk'
      include 'dbnums.blk'

      DIMENSION A(*)
      CHARACTER TYP
      REAL XLN(NPTS,NLNCRV), YLN(NPTS,NLNCRV), ZLN(NPTS,NLNCRV)
      REAL DATA(NWRDS)

      LOGICAL NEED
      CHARACTER T

      CALL DBVIX_BL (TYP, 1, ISID)
      CALL DBVIX_BL (TYP, NVAR, IEID)

      IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN

C      --Read history or global variables

         CALL GETVAR (A, ISID, -999, ISTEP, NWRDS, DATA)

C      --Scan plot variable information and store if a history or global
C      --point from this data is to be plotted

         DO 100 NP = 1, NLNCRV
            CALL DBVTYP_BL (ILVID(1,NP), T, IDUM)
            IF ((T .EQ. 'H') .OR. (T .EQ. 'G')) THEN
               CALL DBVTYP_BL (ILVID(1,NP), T, NE)
               XLN(NPT,NP) = DATA(NE)
               CALL DBVTYP_BL (ILVID(2,NP), T, NE)
               YLN(NPT,NP) = DATA(NE)
               IF (NDIM .GE. 3) THEN
                  CALL DBVTYP_BL (ILVID(3,NP), T, NE)
                  ZLN(NPT,NP) = DATA(NE)
               END IF
            END IF
  100    CONTINUE

      ELSE
         DO 130 ID = ISID, IEID

C         --Determine if variable is needed

            NEED = .FALSE.
            DO 110 NP = 1, NLNCRV
               NEED = NEED .OR. (LOCINT (ID, NDIM, ILVID(1,NP)) .GT. 0)
  110       CONTINUE

            IF (NEED) THEN

C            --Read nodal/element variable

               CALL GETVAR (A, ID, -1, ISTEP, NWRDS, DATA)

C            --Scan plot variable information and store if a node/element
C            --point from this data is to be plotted

               DO 120 NP = 1, NLNCRV
                  NE = ILVNE(NP)
                  IF (ID .EQ. ILVID(1,NP)) XLN(NPT,NP) = DATA(NE)
                  IF (ID .EQ. ILVID(2,NP)) YLN(NPT,NP) = DATA(NE)
                  IF (NDIM .GE. 3) THEN
                     IF (ID .EQ. ILVID(3,NP)) ZLN(NPT,NP) = DATA(NE)
                  END IF
  120          CONTINUE
            END IF

  130    CONTINUE
      END IF

      RETURN
      END
