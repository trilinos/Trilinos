C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LOCNAM (NAME, NAMECO, NAMES,
     &   NENUM, TYPEQV, IDEQV, NERR)
C=======================================================================

C   --*** LOCNAM *** (ALGEBRA) Locate equation variables
C   --Written by Amy Gilkey - revised 02/22/88
C   --
C   --LOCNAM determines the type, location, and node/element specifier (if any)
C   --of an equation variable.  The variable type and database index
C   --are returned.  The specifier is checked for validity.
C   --
C   --Parameters:
C   --   NAME - IN - the variable name
C   --   NAMECO - IN - the coordinate names
C   --   NAMES - IN - the history, global, nodal, and element variable names
C   --   NENUM - IN/OUT - the node/element specifier; changed if incorrect
C   --   TYPEQV - OUT - the variable type (as in /VAR../)
C   --   IDEQV - OUT - the variable database index (as in /VAR../)
C   --   NERR - IN/OUT - the number of errors in the equation, may be set
C   --
C   --Common Variables:
C   --   Uses IXLHS, NAMVAR, TYPVAR of /VAR../
C   --   Uses NDIM, NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/

      include 'ag_namlen.blk'
      include 'ag_var.blk'
      include 'ag_dbnums.blk'

      CHARACTER*(*) NAME
      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      CHARACTER TYPEQV

      CHARACTER*5 STRA, STRB

C   --Save the node/element specifier, if any
      INENUM = NENUM

C   --Search LHS array first; thus, once an input variable is
C   --assigned, its assigned value is referenced

C   --Calculate number of LHS variables
      NUMLHS = MAXVAR - IXLHS + 1

      NV = LOCSTR (NAME, NUMLHS, NAMVAR(IXLHS))
      IF (NV .GT. 0) THEN
         NV = IXLHS + NV - 1
         TYPEQV = TYPVAR(NV)
C      --IDEQV used to determine if referencing a LHS or input variable
         IDEQV = -NV
C      --NENUM slot determined when storage linked
         IF ((TYPEQV .EQ. 'T')
     &      .OR. (TYPEQV .EQ. 'G')) NENUM = -999

      ELSE

C      --Search TIME

         IF (NAME .EQ. 'TIME') THEN
            TYPEQV = 'T'
C         --IDEQV used to determine if referencing a LHS or input variable
            IDEQV = 0
C         --NENUM slot determined when storage linked
            NENUM = -999

         ELSE

C         --Search coordinate names
            NV = LOCSTR (NAME, NDIM, NAMECO)
            IF (NV .GT. 0) THEN
               TYPEQV = 'C'
               IDEQV = NV

            ELSE

C            --Search variable names

               NV = LOCSTR (NAME, NVARGL+NVARNP+NVAREL, NAMES)
               IF (NV .GT. 0) THEN
                  CALL DBVTYP (NV, TYPEQV, ID)
                  IF (TYPEQV .EQ. 'G') THEN
C                  --global variable
                     IDEQV = 1
C                  --NENUM slot determined when storage linked,
C                  --mark with global variable number
                     NENUM = ID
                  ELSE IF (TYPEQV .EQ. 'N') THEN
C                  --Nodal variable
                     IDEQV = ID
                  ELSE IF (TYPEQV .EQ. 'E') THEN
C                  --Element variable
                     IDEQV = ID
                  END IF

               ELSE

C               --Name not found, set up like history to eliminate problems

                  TYPEQV = ' '
                  IDEQV = 0
                  NENUM = 0
                  NERR = NERR + 1
                  CALL PRTERR ('CMDSPEC', '"' // NAME(:LENSTR(NAME))
     &               // '" is not a database or assigned variable')
               END IF
            END IF
         END IF
      END IF

C   --Check validity of node/element specifier

      IF ((TYPEQV .EQ. 'T')
     &   .OR. (TYPEQV .EQ. 'G')) THEN
         IF (INENUM .NE. 0) THEN
            NERR = NERR + 1
            CALL PRTERR ('CMDSPEC', 'Node/element specifier not allowed'
     &         // ' on TIME or global variable '
     &         // NAME(:LENSTR(NAME)))
         END IF

      ELSE IF ((TYPEQV .EQ. 'C') .OR. (TYPEQV .EQ. 'N')) THEN
         IF (NENUM .GT. NUMNP) THEN
            NERR = NERR + 1
            CALL INTSTR (1, 0, NENUM, STRA, LSTRA)
            CALL INTSTR (1, 0, NUMNP, STRB, LSTRB)
            CALL PRTERR ('CMDSPEC', 'Node specifier ' // STRA(:LSTRA)
     &         // ' of ' // NAME(:LENSTR(NAME))
     &         // ' is too big, number of nodes = ' // STRB(:LSTRB))
         END IF

      ELSE IF (TYPEQV .EQ. 'E') THEN
         IF (NENUM .GT. NUMEL) THEN
            NERR = NERR + 1
            CALL INTSTR (1, 0, NENUM, STRA, LSTRA)
            CALL INTSTR (1, 0, NUMNP, STRB, LSTRB)
            CALL PRTERR ('CMDSPEC', 'Element specifier ' // STRA(:LSTRA)
     &         // ' of ' // NAME(:LENSTR(NAME))
     &         // ' is too big, number of elements = ' // STRB(:LSTRB))
         END IF
      END IF

      RETURN
      END
