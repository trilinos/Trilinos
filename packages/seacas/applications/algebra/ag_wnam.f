C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE WNAM (NDBOUT, NDIM, NELBLK, NELBO, VISELB,
     &           NVARGO, NVARNO, NVAREO, NAMECO, BLKTYP,
     &           NAMEGV, NAMENV, NAMEEV, IEVVAR, ISEVOK,
     &           IDELB, NEWID, IEVOK)

C=======================================================================

C   --*** WNAM *** Write the names of the coordinates, element block
C   --             types, and the database result variables. Also,
C   --             write the element block variable truth table
C   --   Modified from RWNAM
C   --*** RWNAM *** (ALGEBRA) Read and write database names
C   --   Written by Amy Gilkey - revised 12/14/87
C   --
C   --Parameters:
C   --   NDBOUT - IN - the output database file
C   --   NDIM   - IN - the number of coordinates per node
C   --   NELBLK - IN - the number of input element blocks
C   --   NELBO  - IN - the number of output element blocks
C   --   VISELB - IN - true iff element block i is to be written
C   --   NVARGO - IN - the number of output global variables
C   --   NVARNO - IN - the number of output nodal variables
C   --   NVAREO - IN - the number of output element variables
C   --   NAMECO - IN - the names of the coordinates
C   --   BLKTYP - IN - the names of the element block types
C   --   NAMEGV - IN - the names of the output global variables
C   --   NAMENV - IN - the names of the output nodal variables
C   --   NAMEEV - IN - the names of the output element variables
C   --   IEVVAR - IN - the ISEVOK index of the output element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --                 variable i of block j exists iff ISEVOK(j,i)
C   --   IDELB  - IN - element block IDs
C   --   IEVOK  - IN -
C   --   NEWID  - IN -
C   --   IOERR  - OUT - input/output error flag
C   --
C   --Common Variables:
C   --   Uses ICOBEG, ICOEND of /DBXVAR/
C   --   Uses NDBIN of /DBASE/
C   --   Uses NDIM, NUMNP, NUMEL, NSTEPS of /DBNUMS/

      include 'exodusII.inc'
      include 'ag_namlen.blk'

      LOGICAL VISELB(NELBLK)
      CHARACTER*(namlen) NAMECO(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(*) NAMEGV(*)
      CHARACTER*(*) NAMENV(*)
      CHARACTER*(*) NAMEEV(*)
      INTEGER IEVVAR(*)
      LOGICAL ISEVOK(NELBLK,NVAREO)
      INTEGER IDELB(*)
      INTEGER IEVOK(NVAREO,NELBO)
      INTEGER NEWID(*)

C     Write the names of the coordinates to the output file
      call expcon(ndbout, nameco, ierr)

C   --Write the variable names

      IF (NVARGO + NVARNO + NVAREO .GT. 0) THEN
        if (nvargo .gt. 0) then
          call expvp(ndbout, 'g', nvargo, ierr)
          call expvan(ndbout, 'g', nvargo, namegv, ierr)
        end if
        if (nvarno .gt. 0) then
          call expvp(ndbout, 'n', nvarno, ierr)
          call expvan(ndbout, 'n', nvarno, namenv, ierr2)
        end if
        if (nvareo .gt. 0) then
          call expvp(ndbout, 'e', nvareo, ierr)
          call expvan(ndbout, 'e', nvareo, nameev, ierr)
        end if
      END IF

      IF ((NVAREO .GT. 0) .AND. (NELBLK .GT. 0)) THEN
         NO = 0
C        Loop from 1 to number of element blocks
         DO 110 IELB = 1, NELBLK
            IF (VISELB(IELB)) THEN
C              Increment counter for writing element block data
               NO = NO + 1
C              New block id index = element block ID
               NEWID(NO) = IDELB(IELB)
C              Loop from 1 to the number of variables to output
               DO 100 I = 1, NVAREO
                  IEV = IEVVAR(I)
                  IF (ISEVOK(IELB,IEV)) THEN
                     IEVOK(I,NO) = 1
                  ELSE
                     IEVOK(I,NO) = 0
                  END IF
 100           CONTINUE
            end if
 110     CONTINUE

C        Write the element variable truth table to the output file
         CALL EXPVTT (NDBOUT, NELBO, NVAREO, IEVOK, IERR)

      END IF

      RETURN
      END

