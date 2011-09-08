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

C $Id: rdxyz.f,v 1.3 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RWXYZ (NTXT, NDB, NDIM, NUMNP, XN, YN, ZN, NAMECO,
     *  NAMLEN,*)
C=======================================================================

C   --*** RDXYZ *** (TXTEXO) Read database coordinates
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDXYZ reads the coordinate array from the text file.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   XN, YN, ZN - OUT - the coordinates
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      REAL XN(*), YN(*), ZN(*)
      CHARACTER*(NAMLEN) NAMECO(*)
      integer idum(3), kval(3)
      real    rdum(3)
      
      character*512 scratch
      CHARACTER*5 STRA

      INP = 0
      READ (NTXT, *, END=110, ERR=110)
      READ (ntxt, '(A)', END=110, ERR=110) scratch
      idcont = 0
      call ffistr (scratch, 3, idcont, nfield, kval, nameco, ival, rval)

      READ (NTXT, *, END=120, ERR=120)
      DO 100 INP = 1, NUMNP
         IF (NDIM .EQ. 2) THEN
            READ (NTXT, *, END=120, ERR=120) XN(INP), YN(INP)
         ELSE IF (NDIM .EQ. 3) THEN
            READ (NTXT, *, END=120, ERR=120) XN(INP), YN(INP), ZN(INP)
         END IF
  100 CONTINUE

      call expcon (ndb, nameco, ierr)
      call expcor (ndb, xn, yn, zn, ierr)

      RETURN

 110  continue
      CALL PRTERR ('FATAL', 'Reading COORDINATE NAMES')
      RETURN 1
      
  120 CONTINUE
      CALL INTSTR (1, 0, INP, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading COORDINATES for node ' // STRA(:LSTRA))
      RETURN 1
      END
