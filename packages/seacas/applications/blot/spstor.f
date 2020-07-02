C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

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

      CALL DBVIX_BL (TYP, 1, ISID)
      CALL DBVIX_BL (TYP, NVAR, IEID)

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
