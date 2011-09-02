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

C $Id: rdeb1.f,v 1.2 2007/10/17 18:47:21 gdsjaar Exp $
C $Log: rdeb1.f,v $
C Revision 1.2  2007/10/17 18:47:21  gdsjaar
C Added copyright notice to all files.
C
C extexo2 is licensed under the BSD license
C
C Revision 1.1.1.1  1991/02/22 19:22:30  gdsjaar
C TxtExo - Convert Text File to EXODUS File
C
c Revision 1.1  1991/02/22  19:22:29  gdsjaar
c Initial revision
c

C=======================================================================
      SUBROUTINE RDEB1 (NTXT, IELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, natrdm, *)
C=======================================================================

C   --*** RDEB1 *** (TXTEXO) Read database element block misc.
C   --   Written by Amy Gilkey - revised 09/30/87
C   --
C   --RDEB1 reads the element block connectivity and attribute information
C   --from the text file.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   IELB - IN - the element block number (for errors)
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - OUT - the element connectivity for this block
C   --   ATRIB - OUT - the attributes for this block
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block misc. information
C   --upon entry; upon exit at end of element block misc. information.

      INTEGER LINK(NUMLNK,*)
      REAL ATRIB(natrdm, *)

      CHARACTER*5 STRA, STRB

      NE = 0
      READ (NTXT, *, END=120, ERR=120)
      DO 100 NE = 1, NUMELB
         READ (NTXT, *, END=120, ERR=120) (LINK(I,NE), I=1,NUMLNK)
  100 CONTINUE

      IF (NUMATR .GT. 0) THEN
         NE = 0
         READ (NTXT, *, END=130, ERR=130)
         DO 110 NE = 1, NUMELB
            READ (NTXT, *, END=130, ERR=130) (ATRIB(I,NE), I=1,NUMATR)
  110    CONTINUE
      END IF

      RETURN

  120 CONTINUE
      CALL INTSTR (1, 0, NE, STRA, LSTRA)
      CALL INTSTR (1, 0, IELB, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading NODES for element ' // STRA(:LSTRA)
     &   // ' for BLOCK ' // STRB(:LSTRB))
      GOTO 140
  130 CONTINUE
      CALL INTSTR (1, 0, NE, STRA, LSTRA)
      CALL INTSTR (1, 0, IELB, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading ATTRIBUTES for element ' // STRA(:LSTRA)
     &   // ' for BLOCK ' // STRB(:LSTRB))
      GOTO 140
  140 CONTINUE
      RETURN 1
      END
