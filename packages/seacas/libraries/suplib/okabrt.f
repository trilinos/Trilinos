C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION OKABRT (ISOK)
C=======================================================================
C$Id: okabrt.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: okabrt.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1990/11/30 09:51:00  gdsjaar
CModified to work on Unicos
C
c Revision 1.1.1.1  90/08/14  16:15:56  gdsjaar
c Testing
c
c Revision 1.1  90/08/14  16:15:54  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:41  gdsjaar
c Initial revision
c

C   --*** OKABRT *** (ETCLIB) Initialize cancel function
C   --   Written by Amy Gilkey - revised 12/21/87
C   --
C   --OKABRT initializes the cancel flag.  It must be called before ISABRT.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag

      LOGICAL ISABRT
      LOGICAL ISOK

      LOGICAL CPUIFC, LDUM

      LOGICAL DOABRT
      SAVE DOABRT

      DATA DOABRT / .FALSE. /

C   --Initialize enable cancel flag
      DOABRT = ISOK

      IF (DOABRT) THEN
C      --Initialize cancel flag
         LDUM = CPUIFC (.TRUE.)
      END IF

      OKABRT = DOABRT

      RETURN

C=======================================================================
      ENTRY ISABRT ()
C=======================================================================
C$Id: okabrt.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: okabrt.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1990/11/30 09:51:00  gdsjaar
CModified to work on Unicos
C
c Revision 1.1.1.1  90/08/14  16:15:56  gdsjaar
c Testing
c
c Revision 1.1  90/08/14  16:15:54  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:41  gdsjaar
c Initial revision
c

C   --*** ISABRT *** (ETCLIB) Check cancel function
C   --   Written by Amy Gilkey - revised 12/17/87
C   --
C   --ISABRT checks the cancel flag.  If it is set, it aborts the current
C   --processing.  In any case, the value of the cancel flag is returned
C   --as the function value.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag

      IF (DOABRT) THEN
C      --Return cancel flag
         ISABRT = CPUIFC (.FALSE.)

         IF (ISABRT) THEN
C         --Send abort message
            WRITE (*, '(1X, A)') '*** Processing aborted ***'
         END IF

      ELSE
         ISABRT = .FALSE.
      END IF

      RETURN
      END
