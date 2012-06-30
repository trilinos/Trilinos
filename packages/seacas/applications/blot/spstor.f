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

C $Log: spstor.f,v $
C Revision 1.3  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1999/03/09 19:41:43  gdsjaar
C Fixed missing parameter definition of MSHBOR in blotII2.f
C
C Cleaned up parameters and common blocks in other routines.
C
C Revision 1.1  1994/04/07 20:15:05  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:33  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SPSTOR (A, ISTEP, TYP, NWRDS, NVAR, NENUM,
     &   PLTVAL, DATA)
C=======================================================================

C   --*** SPSTOR *** (SPLOT) Read and store plot variables from database
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --SPSTOR reads variables from the database and stores any that are
C   --plot variables in the appropriate location.  It reads only one
C   --variable type (nodal or element) per call.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the time step number
C   --   TYP - IN - the variable type ("N"odal or "E"lement)
C   --   NWRDS - IN - the number of words in a data record
C   --   NVAR - IN - the number of records to read
C   --   NENUM - IN - the selected node/element numbers
C   --   PLTVAL - OUT - the plot data array
C   --   DATA - SCRATCH - size = NWRDS
C   --
C   --Common Variables:
C   --   Uses NELBLK, NVAREL of /DBNUMS/
C   --   Uses NNENUM of /SELNE/
C   --   Uses NSPVAR, ISVID of /SPVARS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'selne.blk'
      include 'spvars.blk'

      DIMENSION A(*)
      CHARACTER TYP
      INTEGER NENUM(NNENUM)
      REAL PLTVAL(NNENUM,NSPVAR)
      REAL DATA(NWRDS)

      LOGICAL NEED

      CALL DBVIX (TYP, 1, ISID)
      CALL DBVIX (TYP, NVAR, IEID)

      DO 120 ID = ISID, IEID

C      --Determine if variable is needed

         NEED = (LOCINT (ID, NSPVAR, ISVID) .GT. 0)

         IF (NEED) THEN

C         --Read nodal/element variable

            CALL GETVAR (A, ID, -1, ISTEP, NWRDS, DATA)

C         --Scan plot variable information and store if a node/element
C         --point from this data is to be plotted

            DO 110 N = 1, NSPVAR
               IF (ID .EQ. ISVID(N)) THEN
                  DO 100 I = 1, NNENUM
                     NE = NENUM(I)
                     PLTVAL(I,N) = DATA(NE)
  100             CONTINUE
               END IF
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
