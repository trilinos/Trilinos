C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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

C=======================================================================
      SUBROUTINE RDELB (A, IDELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, KATRIB, *)
C=======================================================================

C   $Id: rdelb.f,v 1.3 1991/01/09 12:59:34 gdsjaar Exp $
C   $Log: rdelb.f,v $
C   Revision 1.3  1991/01/09 12:59:34  gdsjaar
C   Initial conversion from GEN3D to GENSHELL, no BC yet
C
c Revision 1.2  90/10/01  15:41:01  gdsjaar
c Removed MAX(NUMEL,1) from dimension statement.  Note that we now
c assume that NUMEL is at least 1.  This was Non-ANSI usage.
c 
c Revision 1.1.1.1  90/08/20  12:22:40  gdsjaar
c Gen3D Mesh Generation Program
c 
c Revision 1.1  90/08/20  12:22:39  gdsjaar
c Initial revision
c 

C   --*** RDELB *** (GEN3D) Read database element blocks
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --RDELB reads the element block information from the database.
C   --Some dynamic dimensioning is done.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IDELB - OUT - the element block ID for each block
C   --   NUMELB - OUT - the number of elements for each block
C   --   NUMLNK - OUT - the number of nodes per element for each block
C   --   NUMATR - OUT - the number of attributes for each block
C   --   LINK - OUT - the connectivity array (4 nodes per element)
C   --   KATRIB - OUT - the dynamic memory pointer to the attribute array
C   --      (named 'ATRIB', packed)
C   --   * - return statement if end of file or read error
C   --
C   --Common Variables:
C   --   Uses NDBIN of /DBASE/
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --
C   --Database must be positioned at start of element block information
C   --upon entry; upon exit at end of element block information.

      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbnums.blk'

      DIMENSION A(*)
      INTEGER IDELB(NELBLK)
      INTEGER NUMELB(NELBLK)
      INTEGER NUMLNK(NELBLK)
      INTEGER NUMATR(NELBLK)
      INTEGER LINK(2**NDIM,NUMEL)

      CHARACTER*5 STRA

      IEATR = 0
      CALL MDRSRV ('ATRIB', KATRIB, IEATR)
      IEND = 0

      DO 10 IELB = 1, NELBLK

         READ (NDBIN, IOSTAT=IERR, END=30, ERR=30)
     &      IDELB(IELB), NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB)

         ISTART = IEND + 1
         IEND = ISTART + NUMELB(IELB) - 1
         IF (IEND .GT. NUMEL) GOTO 40

         ISATR = IEATR + 1
         IEATR = ISATR + NUMATR(IELB) * NUMELB(IELB) - 1
         CALL MDLONG ('ATRIB', KATRIB, IEATR)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 20

         CALL RELBLK (IELB, NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB),
     &      LINK(1,ISTART), A(KATRIB+ISATR-1), *50)
   10 CONTINUE

   20 CONTINUE
      RETURN

   30 CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading ELEMENT BLOCK SIZING PARAMETERS for block '
     &   // STRA(:LSTRA))
      GOTO 50
   40 CONTINUE
      CALL PRTERR ('FATAL', 'Number of elements in blocks > NUMEL')
      GOTO 50
   50 CONTINUE
      RETURN 1
      END
