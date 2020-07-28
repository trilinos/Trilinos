C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRTDEB (OPTION, NUM)
C=======================================================================

C   --*** PRTDEB *** (ALGEBRA) Print debug information
C   --   Written by Amy Gilkey - revised 12/14/87
C   --
C   --PRTDEB prints out the contents of the common requested.
C   --
C   --Parameters:
C   --   OPTION - IN - the information requested
C   --   NUM - IN - dependent on OPTION
C   --
C   --Common Variables:
C   --   Uses NUMEQN, NUMENT, NAMENT, TYPENT, INXENT, VALENT, VSZENT of /ENT../
C   --   Uses NUMINP, IXLHS, NAMVAR, TYPVAR, IDVAR, ISTVAR of /VAR../

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_ent.blk'
      include 'ag_var.blk'
      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      CHARACTER*(*) OPTION
      CHARACTER TYPE

      IF (OPTION .EQ. 'EQUATION') THEN
         TYPE = 'E'
      ELSE IF (OPTION .EQ. 'VARIABLE') then
         TYPE = 'V'
      ELSE
         WRITE (*, *)
         WRITE (*, 10040) 'DEBUG - Unknown OPTION ', OPTION
         TYPE = ' '
      END IF

      IF (TYPE .EQ. 'E') THEN
         IF (NUM .LE. 0) THEN
            ISTART = 1
            IEND = NUMEQN
         ELSE
            ISTART = NUMEQN
            IEND = NUMEQN
         END IF
         DO 110 NEQN = ISTART, IEND
            WRITE (*, *)
            WRITE (*, 10050, IOSTAT=IDUM) 'Equation ', NEQN
            WRITE (*, 10000)
            DO 100 I = 1, NUMENT(NEQN)
               WRITE (*, 10010, IOSTAT=IDUM) I,
     &            NAMENT(I,NEQN), TYPENT(I,NEQN),
     &            INXENT(I,NEQN), VALENT(I,NEQN),
     &            ITMENT(I,NEQN), IEVENT(I,NEQN), VSZENT(I,NEQN)
10000           FORMAT (4X, ' #', 3X, 'name    ', 3X, 'type',
     &            3X, 'index', 3X, '     value',
     &            3X, 'time', 3X, 'evok', 3X, 'size')
10010           FORMAT (4X, I2, 3X, A, 3X, 3X, A,
     &            3X, I5, 3X, F10.3,
     &            3X, I4, 3X, I4, 3X, 3X, A)
  100       CONTINUE
  110    CONTINUE

      ELSE IF (TYPE .EQ. 'V') THEN
         WRITE (*, *)
         WRITE (*, 10050) 'INPUT variables'
         WRITE (*, 10020)
         DO 120 I = 1, NUMINP
            WRITE (*, 10030, IOSTAT=IDUM) I,
     &         NAMVAR(I), TYPVAR(I), IDVAR(I), (ISTVAR(K,I), K=1,3),
     &         IEVVAR(I)
10020        FORMAT (4X, '  #', 3X, 'name    ', 3X, 'type',
     &         3X, '   id', 3X, '       store   ', 3X, 'evok')
10030        FORMAT (4X, I3, 3X, A, 3X, 3X, A,
     &         3X, I5, 3X, 3I5, 3X, I4)
  120    CONTINUE
         WRITE (*, 10050) 'LHS variables'
         WRITE (*, 10020)
         DO 130 I = IXLHS, MAXVAR
            WRITE (*, 10030, IOSTAT=IDUM) I,
     &         NAMVAR(I), TYPVAR(I), IDVAR(I), (ISTVAR(K,I), K=1,3),
     &         IEVVAR(I)
  130    CONTINUE
      END IF

      IF (TYPE .NE. ' ') WRITE (*, *)

      RETURN
10040  FORMAT (1X, 5A)
10050  FORMAT (1X, A, I5)
      END
