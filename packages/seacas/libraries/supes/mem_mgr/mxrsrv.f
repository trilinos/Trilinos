C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MXRSRV (MYV, NAME1, NEWLEN, NEWLOC, MYLOC, OFFSET,
     *   VOID, LVOID,
     *   NVOIDS, DICT, DPOINT, LDICT, NNAMES, CHRCOL,
     *   DEFER, FILL, FDATA,
     *   LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     This routine finds space to service a non-negative space request.
C     If zero space is requested, a valid pointer of 1 will be
C     generated.

C***********************************************************************

C     MYV      Internal reference array.
               DIMENSION MYV(*)
C     NAME1    Name to be inserted in the dictionary
               CHARACTER*8 NAME1
C     NEWLEN   Length of requested storage
C     NEWLOC   Pointer of new storage (returned)
C     MYLOC    Reference address of internal array
C     OFFSET   Offset between internal array and user's array
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary tables
C     NNAMES   Number of names
               CHARACTER DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     DEFER    Flag for deferred mode.
               LOGICAL DEFER
C     FILL     Flag for data fill.
C     FDATA    Data for fill.
               LOGICAL FILL
C     LASTER   Error return

C***********************************************************************

      LASTER = SUCESS
      MYLEN = NEWLEN

      IF (NEWLEN .EQ. 0) THEN

C        Zero length entry.

         NEWLOC = 1 - OFFSET

      ELSE IF (DEFER) THEN

         CALL MXLOOK (MYLEN, VOID, CHRCOL*LVOID, NVOIDS(1),
     *      VROW, LASTER)

         IF (LASTER .EQ. SUCESS) THEN
            NEWLOC = VOID(VROW,1,1)
         ELSE IF (LASTER .EQ. NOGET) THEN

C           A good void was not found - defer the space request.

            NEWLOC = IXLNUM(NEWLOC)
            MYLEN = - NEWLEN
            LASTER = SUCESS

         END IF

      ELSE

C        Get space.

         CALL MXGET (MYLOC, MYLEN, VOID, LVOID, NVOIDS,
     *      CHRCOL, LASTER, VROW)
         IF (LASTER .NE. SUCESS) RETURN

         NEWLOC = VOID(VROW,1,1)

      END IF

C     Update dictionary.

      CALL MXNSRT (NAME1, NEWLOC, MYLEN, DICT, DPOINT, LDICT,
     *   NNAMES, CHRCOL, LASTER)
      IF (LASTER .EQ. WRTYPE) LASTER = BDNAME
      IF (LASTER .NE. SUCESS) RETURN

      IF (MYLEN .GT. 0) THEN

C        Data fill pattern.

         IF (FILL) THEN
            DO 100 I = VOID(VROW,1,1), VOID(VROW,1,1)+MYLEN-1-7,8
               MYV(I+0) = FDATA
               MYV(I+1) = FDATA
               MYV(I+2) = FDATA
               MYV(I+3) = FDATA
               MYV(I+4) = FDATA
               MYV(I+5) = FDATA
               MYV(I+6) = FDATA
               MYV(I+7) = FDATA
  100       CONTINUE
            do 110 J = I, VOID(VROW,1,1)+MYLEN-1
              MYV(J) = FDATA
 110        continue
         END IF

C        Update void table.

         VOID(VROW,1,1) = VOID(VROW,1,1) + MYLEN
         VOID(VROW,1,2) = VOID(VROW,1,2) - MYLEN
         CALL VTABLE (1, 0, VOID, LVOID, NVOIDS(1), CHRCOL, LASTER)
         NEWLOC = NEWLOC + OFFSET

      ELSE IF (MYLEN .LT. 0) THEN

         NEWLOC = OFFSET - MYLOC

      END IF

      RETURN
      END
