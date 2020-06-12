C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CKWHOL (WHOTIM, *)
C=======================================================================
C$Id: ckwhol.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: ckwhol.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:13  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:12  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:05  gdsjaar
c Initial revision
c

C   --*** CKWHOL *** (ETCLIB) Check for whole time step
C   --   Written by Amy Gilkey - revised 12/23/87
C   --
C   --CKWHOL prints an error message if the current time step is not
C   --a whole time step.
C   --
C   --Parameters:
C   --   WHOTIM - IN - true iff whole time step
C   --   * - return statement if error

      LOGICAL WHOTIM

      IF (.NOT. WHOTIM) THEN
         CALL PRTERR ('CMDERR',
     &      'Only history variables are present on this time step')
         RETURN 1
      END IF

      RETURN
      END
