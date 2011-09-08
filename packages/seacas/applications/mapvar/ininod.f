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

      SUBROUTINE ININOD(SOLNB,IVAR,TIMES,ISTP,IDBLK,NDLSTB,XB,YB,ZB,
     &                  SN)
C
C *********************************************************************
C
C ININOD initializes nodal variable values based on TIME, ELEMENT
C BLOCK, VARIABLE NAME, COORDINATE, etc. By default, nodal variable 
C values are set  to zero. It is intended that the user rewrite this 
C subroutine to provide values that are appropriate to the problem 
C being solved. This is the preferred method to handle nodal variable 
C assignment for recipient mesh nodes that lie outside the boundary 
C of the donor mesh.
C
C
C Called by INTRPN, SINTPN
C
C *********************************************************************
C
C SOLNB   REAL  Array of nodal variable values 
C               (1:nodesb,1:nvarnp,1:ntimes)
C IVAR    INT   Position of variable in SOLNB array
C TIMES   REAL  Array of times (1:ntimes)
C ISTP    INT   Position in TIMES array, time step being processed
C IDBLK   INT   The element block I. D.
C NDLSTB  INT   Array of nodes that belong to this element block
C               (1:numndb)
C XB      REAL  Array of X-coordinates (1:nodesb)
C YB      REAL  Array of Y-coordinates (1:nodesb)
C ZB      REAL  Array of Z-coordinates (1:nodesb)
C SN      REAL  Dummy array to store values from MESH-B
C
C *********************************************************************
C
      include 'exodusII.inc'
C
      include 'aexds1.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
C
      DIMENSION SOLNB(NODESB,NVARNP), TIMES(*), NDLSTB(*) 
      DIMENSION XB(*), YB(*), ZB(*)
      DIMENSION SN(*)
C
C *********************************************************************
C
C Code to help you find some potentially useful stuff
C The actual time (real number) 
C     TIME = TIMES(ISTP)
C
C The pointer into VARNAM to get the variable name being processed
C     INAM = IVAR + NVARGP + NVAREL
C
C The name of the variable (character) being processed
C     NAME = NAMVAR(INAM)
C
C INOD = NDLSTB(I)
C
C Set value of node-NDLSTB(I); variable-IVAR to 0. by default.
C User to replace this with whatever code he wishes.
C
      CALL EXGNV(NTP3EX,ISTP,IVAR,NODESB,SN,IERR)
      DO 10 I = 1, NUMNDB
        INODE = NDLSTB(I)
        SOLNB(INODE,IVAR) = SN(INODE)
 10   CONTINUE
C
      RETURN
      END
