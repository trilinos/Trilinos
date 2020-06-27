C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: addv.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C $Log: addv.f,v $
C Revision 1.2  2009/03/25 12:46:01  gdsjaar
C Add copyright and license notice to all files.
C
C Revision 1.1  1993/11/18 21:32:06  gdsjaar
C Added scilib routines saxpy.f scopy.f sdot.f snrm2.f subv.f
C Added utility routines addv.f subv.f
C
C-----------------------------------------------------------------------
      SUBROUTINE ADDV( N,A,B,C )
C
C***********************************************************************
C
C     DESCRIPTION: This routine adds two vectors
C
C     FORMAL PARAMETERS:
C        N        INTEGER   Number of entries in A, B
C        A        REAL      First vector
C        B        REAL      Vector to be added
C        C        REAL      Vector with the result
C
C***********************************************************************
C
      DIMENSION A(N),B(N),C(N)
C
      DO 100 I = 1,N
        C(I) = A(I) + B(I)
  100 CONTINUE
C
      RETURN
      END
