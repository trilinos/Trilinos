C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBIVTT (NDB, ISEVOK, ITMP, NELBLK, NVAREL)
C=======================================================================
C   --*** DBIVTT *** Read element variable truth table
C   --   Modified for ExodusII format 8/26/95
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   ISEVOK - OUT - the dynamic memory array of the element block variable
C   --                  truth table;variable i,block j exists iff ISEVOK(j,i)
C   --   ITMP   - OUT - temporary scratch array
C   --   NELBLK - IN  - the number of element blocks
C   --   NVAREL - IN  - the number of element variables;

      INTEGER NDB
      INTEGER NELBLK, NVAREL
      LOGICAL ISEVOK(NELBLK, NVAREL)
      INTEGER ITMP(NVAREL, NELBLK)

C     Read the element block variable truth table
C     call exgvtt(fileid, num_elem_blks, num_elem_var,
C                 isevok(num_elem_var, num_elem_blks, errorid)
C       isevok - num_elem_var cycles faster
      if (nvarel .gt. 0) then
         CALL EXGVTT(NDB, NELBLK, NVAREL, ITMP, IERR)
         DO 110 I = 1, NVAREL
            DO 100 IELB = 1, NELBLK
              ISEVOK(IELB,I) = (ITMP(I,IELB) .NE. 0)
 100        CONTINUE
 110     CONTINUE
      end if

      RETURN
      END

