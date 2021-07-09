C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBLIST (TYPE, A, NAMECO, BLKTYP, NAMES,
     &                   TIMES, IDELB, QAREC, INFREC, MERR)
C=======================================================================
c      SUBROUTINE DBLIST (TYPE, A,
c     &   NAMECO, BLKTYP, NAMES, TIMES, WHOTIM, IDELB, QAREC, INFREC)

C   --*** DBLIST *** (ALGEBRA) Display database information
C   --   Written by Amy Gilkey - revised 01/07/88
C   --
C   --DBLIST prints out the requested information.
C   --
C   --Parameters:
C   --   TYPE   - IN - the type of LIST requested (none, VARS)
C   --   A      - IN - the dynamic memory base array
C   --   NAMECO - IN - the coordinate names
C   --   BLKTYP - IN - the element block names
C   --   NAMES  - IN - the global, nodal, and element variable names
C   --   TIMES  - IN - the database time steps
C   --   IDELB  - IN - the element block IDs
C   --   QAREC  - IN - the QA records containing:
C   --           (1) - the analysis code name
C   --           (2) - the analysis code QA descriptor
C   --           (3) - the analysis date
C   --           (4) - the analysis time
C   --   INFREC - IN - the information records
C   --
C   --Common Variables:
C   --   Uses NDBIN of /DBASE/
C   --   Uses TITLE of /DBTITL/
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK, NVARNP, NVAREL, NVARGL
C   --      of /DBNUMS/
C   --   Uses NQAREC, NINFO of /DBNUMQ/

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_dbase.blk'
      include 'ag_dbtitl.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbnumg.blk'
      include 'ag_dbnumq.blk'

      CHARACTER*(*) TYPE
      DIMENSION A(*)
      CHARACTER*(namlen) NAMECO(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(namlen) NAMES(*)
      REAL TIMES(*)
      INTEGER IDELB(*)
      CHARACTER*(MXSTLN) QAREC(4,*)
      CHARACTER*(MXLNLN) INFREC(*)

      CHARACTER*(mxstln) SHOTYP
      CHARACTER*(MXLNLN) STRING
      LOGICAL LDUM
      CHARACTER CDUM

      CHARACTER*(mxstln) SHOTBL(8)
      INTEGER MERR
      SAVE SHOTBL
C      --SHOTBL - the DBLIST type table

      MERR = 0

C   --DBLIST type table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA SHOTBL /
     1  'VARS                            ',
     *  'BLOCKS                          ',
     *  'MATERIAL                        ',
     *  'QA                              ',
     *  'NAMES                           ',
     2  'STEPS                           ',
     *  'TIMES                           ',
     3  '                                ' /

      CALL ABRSTR (SHOTYP, TYPE, SHOTBL)
      IF (SHOTYP .EQ. ' ') SHOTYP = TYPE

      IF (SHOTYP .EQ. 'VARS') THEN
          CALL DBPINI ('TISV', NDBIN, TITLE,
     &         NDIM, NUMNP, NUMEL, NELBLK,
     &         NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL,
     &         LESSDF, NVARGL, NVARNP, NVAREL, ' ')

      ELSE IF ((SHOTYP .EQ. 'BLOCKS')
     &    .OR. (SHOTYP .EQ. 'MATERIAL')) THEN
          CALL MDFIND ('NUMELB', KNELB, IDUM)
          CALL MDFIND ('NUMLNK', KNLNK, IDUM)
          CALL MDFIND ('NUMATR', KNATR, IDUM)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
             CALL MEMERR
             MERR = 1
          END IF
          CALL DBPELB ('N', NELBLK, IDELB,
     &         A(KNELB), A(KNLNK), A(KNATR),
     &         BLKTYP, IDUM, CDUM, LDUM, IDUM)

      ELSE IF (SHOTYP .EQ. 'QA') THEN
          CALL DBPQA ('*', NQAREC, QAREC, NINFO, INFREC)

      ELSE IF (SHOTYP .EQ. 'NAMES') THEN
          WRITE (*, *)
          WRITE (STRING, 10000) (NAMECO(I), I=1,NDIM)
10000     FORMAT (6(' ',A8))
          CALL SQZSTR (STRING, LSTR)
          WRITE (*, 10010) 'Coordinate names: ', STRING(:LSTR)
10010     FORMAT (1X, 10A)

          CALL DBVIX ('G', 1, IGV)
          CALL DBVIX ('N', 1, INV)
          CALL DBVIX ('E', 1, IEV)
          CALL DBPNAM ('*', NVARGL, NVARNP, NVAREL,
     &         NAMES(IGV), NAMES(INV),NAMES(IEV))

      ELSE IF (SHOTYP .EQ. 'STEPS') THEN
          CALL DBPTIM ('NM', NSTEPS, TIMES)

      ELSE IF (SHOTYP .EQ. 'TIMES') THEN
          CALL DBPTIM ('NT', NSTEPS, TIMES)

      ELSE
          CALL SHOCMD ('LIST Options:', SHOTBL)
      END IF

      RETURN
      END
