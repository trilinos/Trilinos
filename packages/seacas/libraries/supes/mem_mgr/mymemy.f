C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYMEMY( MEMREQ, LOCBLK, MEMRTN, MAXSIZ )
      SAVE NUSED
      DATA NUSED /0/

************************************************************************

C     FORTRAN EXTENSION LIBRARY - ANSI FORTRAN - USER INTERFACE ROUTINE

C     DESCRIPTION:
C     This routine requests the operating system to allocate or release
C     numeric storage. A positive MEMREQ indicates a request for memory,
C     while a negative MEMREQ indicates a release. All locations and
C     sizes are measured in numeric storage units.

C     In memory request mode, MEMRTN .LT. MEMREQ indicates an error.

C     In memory release mode, MEMRTN .LE. -MEMREQ. Furthermore, memory
C     must be released from the top down, i.e., LOCBLK must not change.

C     This version actually allocates storage from a static pool, whose
C     size is defined by the parameter MAXSIZ. If system dependent
C     support for the function IXLNUM is not implemented, the PARAMETER
C     and COMMON statements above must be duplicated in the caller.

C     FORMAL PARAMETERS:
C     MEMREQ    INTEGER         Number of numeric units
C     LOCBLK    INTEGER         Location of memory block
C     MEMRTN    INTEGER         Size of memory block at routine completion
C     MAXSIZ    INTEGER         Size of character memory - dimension in
C                               MDINIT.

C     SAVED VARIABLES:
C     NUSED     INTEGER         Number of units dynamically allocated

************************************************************************

      IF ( MEMREQ .GE. 0 ) THEN

C Allocate storage -
         LOCBLK = 1 + NUSED
         MEMRTN = MIN( MAXSIZ-NUSED , MEMREQ )
         NUSED = NUSED + MEMRTN
      ELSE

         MEMRTN = -MEMREQ
      END IF

      RETURN
      END
