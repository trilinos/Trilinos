C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: nwhsel.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:06:19  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:20  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      INTEGER FUNCTION NWHSEL (NPTIMS, IPTIMS, WHOTIM)
C=======================================================================

C   --*** NWHSEL *** (TIMSEL) Return number of whole times selected
C   --   Written by Amy Gilkey - revised 11/11/87
C   --
C   --NWHSEL returns the number of whole times selected.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of selected times
C   --   IPTIMS - IN - the selected time steps
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step

      INTEGER NPTIMS
      INTEGER IPTIMS(*)
      LOGICAL WHOTIM(*)

      NWHSEL = 0
      DO 100 I = 1, NPTIMS
         IF (WHOTIM(IPTIMS(I))) NWHSEL = NWHSEL + 1
  100 CONTINUE

      RETURN
      END
