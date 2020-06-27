C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: ldrota.f,v 1.1 1993/07/16 16:46:34 gdsjaar Exp $
C $Log: ldrota.f,v $
C Revision 1.1  1993/07/16 16:46:34  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE LDROTA(AXIS,COSANG,SINANG,MAT)
      REAL MAT(4,4)
      CHARACTER*(*) AXIS
      CHARACTER*1 TAXIS

      CALL MXIDEN(4,MAT)
      TAXIS = AXIS
      CALL CHRUP(TAXIS,TAXIS)
      IF (TAXIS.EQ.'X') THEN
         MAT(2,2) = COSANG
         MAT(2,3) = SINANG
         MAT(3,2) = -SINANG
         MAT(3,3) = COSANG

      ELSE IF (TAXIS.EQ.'Y') THEN
         MAT(1,1) = COSANG
         MAT(1,3) = -SINANG
         MAT(3,1) = SINANG
         MAT(3,3) = COSANG

      ELSE IF (TAXIS.EQ.'Z') THEN
         MAT(1,1) = COSANG
         MAT(1,2) = SINANG
         MAT(2,1) = -SINANG
         MAT(2,2) = COSANG
      END IF

      RETURN

      END
