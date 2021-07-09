C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRNODE (OPTION, NOUT, NUMNP, LISNP,
     &   NVARNP, LISNV, NAMENV, VARNP, MAPND)
C=======================================================================

C   --*** PRNODE *** (BLOT) Display current database nodal variables
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --PRNODE displays the nodal data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMNP - IN - the number of nodes
C   --   LISNP - IN - the indices of the selected nodes
C   --   NVARNP - IN - the number of nodal variables
C   --   LISNV - IN - the indices of the selected nodal variables
C   --   NAMENV - IN - the names of the nodal variables
C   --   VARNP - IN - the selected nodal variables for the time step

      CHARACTER*(*) OPTION
      INTEGER LISNP(0:*)
      INTEGER LISNV(0:*)
      CHARACTER*(*) NAMENV(*)
      REAL VARNP(NUMNP,*)
      INTEGER MAPND(*)

      LOGICAL ISABRT

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      do i=1, lisnv(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
          WRITE (NOUT, 10010) irow, icol, NAMENV(LISNV(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMENV(LISNV(I))
        END IF
      end do

      DO 100 IX = 1, LISNP(0)
         IF (ISABRT ()) RETURN
         INP = LISNP(IX)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020, IOSTAT=IDUM)
     &       MAPND(INP), (VARNP(INP,I), I=1,LISNV(0))
         ELSE
            WRITE (*, 10020, IOSTAT=IDUM)
     &       MAPND(INP), (VARNP(INP,I), I=1,LISNV(0))
         END IF
  100 CONTINUE

      RETURN

10000 FORMAT (/, 1X, 'NODAL TIME STEP VARIABLES')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
10020 FORMAT (1X, 'Node', I9, 5 (1X, 1PE13.6), :, /,
     &   (1X, 13X, 5 (1X, 1PE13.6)))
      END
