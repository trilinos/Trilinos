C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRGLOB (OPTION, NOUT, NVARGL, LISGV, NAMEGV, VARGL)
C=======================================================================

C   --*** PRGLOB *** (EXPLORE) Display current database global variables
C   --
C   --PRGLOB displays the global data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARGL - IN - the number of global variables
C   --   LISGV - IN - the indices of the selected global variables
C   --   NAMEGV - IN - the names of the global variables
C   --   VARGL - IN - the global variables for the time step

      include 'exodusII.inc'
      CHARACTER*(*) OPTION
      INTEGER LISGV(0:*)
      CHARACTER*(*) NAMEGV(*)
      REAL VARGL(*)
      INTEGER GETPRC, PRTLEN
      CHARACTER*128 FMT1, FMT

      PRTLEN = GETPRC() + 7
      WRITE(FMT1,20) PRTLEN, PRTLEN-7
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMT, 30) FMT1(:LFMT)

      LNAM = 0
      do i=1, lisgv(0)
        L = lenstr(namegv(lisgv(i)))
        lnam = max(l, lnam)
      end do

      IF (NOUT .GT. 0) THEN
        WRITE (NOUT, 10000)
      ELSE
        WRITE (*, 10000)
      END IF

      do 100 i=1, lisgv(0)
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, FMT) NAMEGV(LISGV(I))(:LNAM), VARGL(LISGV(I))
        ELSE
           WRITE (*, FMT)    NAMEGV(LISGV(I))(:LNAM), VARGL(LISGV(I))
        END IF
 100  continue

      RETURN

 20   FORMAT('1PE',I2.2,'.',I2.2)
 30   FORMAT ('(1X, A, '' = '',', A,')')

10000  FORMAT (/, 1X, 'Global Time Step Variables')
      END
