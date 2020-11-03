C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PROCHS (A, IA, NELBLK, IDELB, IDSCR, NELB, NLNK, NATR,
     &           KLINKE, KATRIB, NAMELB, LPTR, ISHEX, HEXID, *)
C=======================================================================
C   --PROCHS - PROCess HexShell elements
C   --To visualize a HEXSHELL element the 12-noded HEXSHELL needs to
C   --be separated into an 8-noded HEX and 4-node SHELL.  A New element
C   --block will be created for every element block that contains a
C   --HEXSHELL element.
C   --An error message is displayed if the end of file is read.
C   --Parameters:
C   --   IA     - I/O - dynamic memory array for integer values
C   --   A      - I/O - dynamic memory array for real values
C   --   C      - I/0 - dynamic memory array for character values
C   --   NELBLK - IN  - number of element blocks
C   --   IDELB  - I/O - array: element block ids
C   --   NELB   - I/O - array: number of elements per element block
C   --   NLNK   - I/O - array: number of nodes per element
C   --   NATR   - I/O - array: number of attributes
C   --   LINK   - I/O - array: link array
C   --   ATRIB  - I/O - array: attribute array
C   --   LPTR   - IN  - array: link pointer array
C   --   ISHEX  - IN  - number of element blocks with hex shells
C   --   KHEXID - IN  - array storing HEXSHELL element block ids

      DIMENSION A(*),IA(*)
      INTEGER NELBLK
      CHARACTER*(*) NAMELB(*)
      INTEGER IDELB(*), IDSCR(*), NELB(*), NLNK(*), NATR(*)
      INTEGER KLINKE, KATRIB
      INTEGER LPTR(*)
      INTEGER ISHEX, HEXID(*)

      IHEX = 1
 100  CONTINUE
      IF (IHEX .GT. ISHEX) GO TO 200
C     HEXSHELL element block ID - will remain HEX element block id
      IDH = HEXID(IHEX)
C     New element block ID for SHELL
      IDS = IDH + 5000
C     Sort the element block id array
      DO 10 J = NELBLK+ISHEX, 1, -1
         IF (IDELB(J) .EQ. IDH) THEN
C           Element block id for shell element block
            IDELB(J+1) = IDS
C           Number of elements in shell element block
            NELB(J+1)  = NELB(J)
C           Number of nodes per element for hex element block
            NLNK(J)    = 8
C           Number of nodes per element for shell element block
            NLNK(J+1)  = 4
C           Number of attributes in hex element block
            NATR(J)    = 0
C           Number of attributes in shell element block
            NATR(J+1)  = 1
C           Named type for hex element block
            NAMELB(J)  = 'HEX'
C           Named type for shell element block
            NAMELB(J+1)= 'SHELL'
C           Process link array
C           Pointer to LINK array for element block J
            CALL GETPTR (IDH,NELBLK,IDSCR,LPTR,IPTR)
C           Size of link array for HEXSHELL element block
            NELEM   = NELB(J)
C           Reserve space for link array - HEX element block
            CALL MDRSRV('HSCR', KHSCR, 8*NELEM)
C           Reserve space for link array - SHELL element block
            CALL MDRSRV('SSCR', KSSCR, 4*NELEM)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) RETURN 1
            CALL MODLNK(NELEM, IA(KHSCR), IA(KSSCR), IA(IPTR))
            CALL MDDEL('HSCR')
            CALL MDDEL('SSCR')
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) RETURN 1
            GOTO 150
         ELSE
            IDELB(J)  = IDELB(J-1)
            NELB(J)   = NELB(J-1)
            NLNK(J)   = NLNK(J-1)
            NAMELB(J) = NAMELB(J-1)
         ENDIF
  10  CONTINUE
 150  CONTINUE
      IHEX = IHEX + 1
      GO TO 100
 200  CONTINUE

      RETURN
      END

      SUBROUTINE GETPTR(ID, NELBLK, IDS, LPTR, IPTR)
C     ID     - IN  - element block id of link array
C     NELBLK - IN  - number of element blocks
C     IDS    - IN  - element blocks id of original elements
C     LPTR   - IN  - 1st pointers in each element block for LINKE array
C     IPTR   - OUT - pointer into link array

      INTEGER ID, NELBLK, IDS(*), LPTR(*)

C     Search for element block id ID in IDS
      DO 400 I = 1, NELBLK
         IF (ID .EQ. IDS(I)) THEN
            IPTR = LPTR(I)
            GOTO 450
         ENDIF
 400  CONTINUE
 450  CONTINUE

      RETURN
      END

      SUBROUTINE MODLNK(NELEM, HEX, SHELL, LINKE)
C     The LINKE array contains the connectivity for a 12 node HEXSHELL element.
C     This subroutine will create a connectivity array for a HEX element block
C     and a connectivity array for a SHELL element block.  To create the
C     HEX connectivity array the first 8 of 12 nodes from every element in
C     the HEXSHELL element block connectivity are copied into the HEX(*)
C     array (e.g. if the HEXSHELL element block has 2 elements, LINKE(1-8)
C     and LINKE(13-20) are copied into the HEX array).  The SHELL connectivity
C     is created by copying nodes 9 through 12 from every element in the
C     HEXSHELL element block connectivity into the SHELL connectivity
C     (e.g. For 2 HEXSHELL elements, LINKE(9-12) and LINKE(21-24) will be
C     copied into the SHELL connectivity). After the HEX and SHELL
C     connectivity arrays have been created, they are copied back into
C     the LINKE array.  First HEX and then SHELL ares copied back into LINKE.
C     C. Forsythe 7/31/97

      INTEGER NELEM
      INTEGER HEX(*), SHELL(*), LINKE(*)
      INTEGER IH, IS

      IH = 8*NELEM
      IS = 4*NELEM
C     Create HEX connectivity array - copy the first 8 of 12 nodes from
C     every HEXSHELL element into the HEX connectivity array.
      DO 300 I = 1, IH
        INC = (I-1)/8
        J   = I + INC*4
        HEX(I) = LINKE(J)
 300  CONTINUE

C     Create SHELL connectivity array - copy nodes 9-12 from every HEXSHELL
C     element into the SHELL connectivity array.
      DO 310 I = 1, IS
         INC = (I-1)/4 + 1
         J   = I + INC*8
         SHELL(I) = LINKE(J)
 310  CONTINUE

C     Copy the HEX and then the SHELL link arrays back into the LINKE array
      DO 320 I = 1, IH
         LINKE(I) = HEX(I)
 320  CONTINUE
      DO 330 I = 1, IS
         J = IH + I
         LINKE(J) = SHELL(I)
 330  CONTINUE

      RETURN
      END

