C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: skprnt.f,v 1.1 1990/11/30 11:15:48 gdsjaar Exp $
C $Log: skprnt.f,v $
C Revision 1.1  1990/11/30 11:15:48  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SKPRNT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SKPRNT (LUNIT, LEN, RSTACK, NDIM)
C***********************************************************************
C
C  SUBROUTINE SKPRNT = STACK PRINTING ROUTINE
C
C***********************************************************************
C
C
      REAL RSTACK (NDIM)
C
      WRITE (LUNIT, '(I8,G12.5)') (I, RSTACK(I + 2), I = LEN, 1, -1)
C
      RETURN
C
      END
