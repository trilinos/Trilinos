C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: sqzixv.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:15:26  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:39  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SQZIXV (NNE, ISEGEL, VALIN, VALOUT)
C=======================================================================

C   --*** SQZIXV *** (SPLOT) Compress values using indices
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --SQZIXV compresses a list of values using a list of indices.
C   --
C   --Parameters:
C   --   NNE - IN - the number of values to be output
C   --   ISEGEL - IN - the indices of the values to be output
C   --   VALIN - IN - the input values
C   --   VALOUT - OUT - the output (compressed) values

      INTEGER ISEGEL(*)
      REAL VALIN(*)
      REAL VALOUT(*)

      DO 100 I = 1, NNE
         VALOUT(I) = VALIN(ISEGEL(I))
  100 CONTINUE

      RETURN
      END
