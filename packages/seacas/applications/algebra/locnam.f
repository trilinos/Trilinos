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

      include 'namlen.blk'
      include 'var.blk'
      include 'dbnums.blk'

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
