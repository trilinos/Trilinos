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
      SUBROUTINE WRELB (A, BLKTYP, 
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB,
     &   IXEL, INCEL, NREL, IELCOL, IXNP, NRNP)
C=======================================================================

C   $Id: wrelb.f,v 1.2 1991/01/09 12:59:45 gdsjaar Exp $
C   $Log: wrelb.f,v $
C   Revision 1.2  1991/01/09 12:59:45  gdsjaar
C   Initial conversion from GEN3D to GENSHELL, no BC yet
C
c Revision 1.1.1.1  90/08/20  12:23:25  gdsjaar
c Gen3D Mesh Generation Program
c 
c Revision 1.1  90/08/20  12:23:23  gdsjaar
c Initial revision
c 

C   --*** WRELB *** (GEN3D) Write 3D element blocks
C   --   Written by Amy Gilkey - revised 09/02/87
C   --
C   --WRELB calculates and writes the element block information for the
C   --3D database.
C   --Some dynamic dimensioning is done.
C   --
C   --Parameters:
C   --   A - IN/OUT - the dynamic memory base array
C   --   BLKTYP - IN - the element block type
C   --   IDELB - IN - the id for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   LINK - IN - the 2D connectivity array (always 4 nodes)
C   --   ATRIB - IN - the 2D attribute array (packed)
C   --   IXEL - IN - the new index for each element
C   --   INCEL - IN - the increment for each element, needed for blocks
C   --      that become multiple blocks
C   --   NREL - IN - the number of new elements generated for each element
C   --   IELCOL - IN - the row number for each element, 0 if not needed
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses NUMEL3 of /DBNUM3/
C   --   Uses NEREPL of /PARAMS/
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'
      INCLUDE 'gs_xyzmir.blk'

      DIMENSION A(*)
      CHARACTER BLKTYP(NELBLK)
      INTEGER IDELB(NELBLK)
      INTEGER NUMELB(NELBLK)
      INTEGER NUMLNK(NELBLK)
      INTEGER NUMATR(NELBLK)
      INTEGER LINK(4,*)
      REAL ATRIB(*)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)
      INTEGER IXNP(*), NRNP(*)

      MAXID = 0
      MAXEL = 0
      MAXATR = 0
      DO 10 IELB = 1, NELBLK
         MAXID = MAX (MAXID, IDELB(IELB))
         MAXEL = MAX (MAXEL, NUMELB(IELB))
         MAXATR = MAX (MAXATR, NUMATR(IELB))
   10 CONTINUE
      CALL MDRSRV ('ATRIB3', KATR3, MAXATR * MAXEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 50

      IEL = 1
      IATR = 1

      DO 40 IELB = 1, NELBLK

C      --Block id and start/end for blocks that become multiple blocks

         NRSTRT = 1
         IDEBL3 = IDELB(IELB)
         IBLK = 1

            NREND = NEREPL
            NBLK = 1
         IF (NBLK .GT. 1) THEN
            LSTID = MAXID
            MAXID = MAXID + NBLK-1
         END IF

   20    CONTINUE

C      --Number of elements in the block

         NELB3 = NUMELB(IELB) * (NREND - NRSTRT + 1)

C      --Number of nodes per element

         NUMLN3 = NUMLNK(IELB) 

C      --Connectivity and attributes

         CALL NEWEL1 (BLKTYP(IELB),
     &        NUMELB(IELB), NELB3, NRSTRT, NREND,
     &        NUMLNK(IELB), NUMLN3, NUMATR(IELB), MAX(NUMATR(IELB),1),
     $        NUMNP, NUMNP3, LINK(1,IEL), A(KLINK3), ATRIB(IATR),
     $        A(KATR3), IXEL(IEL), INCEL(IEL), NREL(IEL), IELCOL(IEL),
     $        IXNP, NRNP)

C      --Fixup connectivity if mirrored

         IF (XMIRR * YMIRR * ZMIRR .LT. 0.0) THEN
            CALL DBMIRR (IELB, IELB, IDEBL3, NELB3, NUMLN3, A(KLINK3))
         END IF

C      --Write 3D

         CALL DBOELB (NDBOUT, IELB, IELB,
     &      IDEBL3, NELB3, NUMLN3, NUMATR(IELB),
     &      A(KLINK3), A(KATR3))

         IF (IBLK .LT. NBLK) THEN
            IBLK = IBLK + 1
            LSTID = LSTID + 1
            IDEBL3 = LSTID
            NRSTRT = NREND + 1
            GOTO 20
         END IF

         IEL = IEL + NUMELB(IELB)
         IATR = IATR + NUMATR(IELB) * NUMELB(IELB)
   40 CONTINUE

   50 CONTINUE
      CALL MDDEL ('LINK3')
      CALL MDDEL ('ATRIB3')

      RETURN
      END
