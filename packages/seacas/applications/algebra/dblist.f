C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

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

      include 'params.blk'
      include 'namlen.blk'
      include 'dbase.blk'
      include 'dbtitl.blk'
      include 'dbnums.blk'
      include 'dbnumg.blk'
      include 'dbnumq.blk'

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
     &         NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &         LESSDF, NVARGL, NVARNP, NVAREL)

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
