C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C $Id: rdelb.f,v 1.4 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RDELB (NTXT, IELB, IDELB, NUMELB, NUMLNK, NUMATR,
     &   NAMELB, A, KLINK, KATRIB, *)
C=======================================================================

C   --*** RDELB *** (TXTEXO) Read database element block
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDELB reads the element block information from the text file.
C   --Some dynamic dimensioning is done.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   IELB - IN - the element block number
C   --   IDELB - IN - the id for this block
C   --   NUMELB - IN - the number of elements in this block
C   --   NUMLNK - IN - the number of nodes per element in this block
C   --   NUMATR - IN - the number of attributes in this block
C   --   A - IN - the dynamic memory base array
C   --   KLINK - IN - pointer to the element connectivity for this block
C   --   KATRIB - IN - pointer to the attributes for this block
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block information
C   --upon entry; upon exit at end of element block information.

      include 'exodusII.inc'
      DIMENSION A(*)
      CHARACTER*(MXSTLN) NAMELB
      CHARACTER*5 STRA

      NAMELB = ' '
      READ (NTXT, *, END=110, ERR=110)
      READ (NTXT, 130, END=110, ERR=110) IDELB, NUMELB, NAMELB
C ... Strip everything in namelb from first space to end
      IEX = index(namelb, " ")
      if (iex .gt. 0) then
        namelb(iex:) = " "
      end if

      READ (NTXT, *, END=110, ERR=110) NUMLNK, NUMATR

      IECON = NUMLNK * NUMELB
      CALL MDLONG ('LINK', KLINK, IECON)
      IEATR = NUMATR * NUMELB
      CALL MDLONG ('ATRIB', KATRIB, IEATR)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      CALL RDEB1 (NTXT, IELB, NUMELB, NUMLNK, NUMATR,
     &   A(KLINK), A(KATRIB), max(1,numatr), *120)

  100 CONTINUE
      RETURN

  110 CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading ELEMENT BLOCK SIZING PARAMETERS for block '
     &   // STRA(:LSTRA))
  120 CONTINUE
 130  format (2I10,6X,A)
      RETURN 1
      END
