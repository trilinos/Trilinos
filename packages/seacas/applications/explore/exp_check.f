C    Copyright(C) 2018 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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
      subroutine check(a, ia, exodus, idelb, ebtype,
     *  numelb, isevok, numlnk,
     *  numatr, link, atrib, atname, mapnd, dbmapnd, mapel, dbmapel,
     *  idnps, nnnps, ixnnps, ltnnps, facnps, idess, neess, nness,
     *  ixeess, ixness, lteess, ltsess, facess, vargl, varnp, varel)
      implicit none

      include 'exodusII.inc'
      INCLUDE 'exp_dbnums.blk'

      REAL A(*)
      INTEGER IA(*)
      LOGICAL EXODUS
      INTEGER IDELB(*), NUMELB(*)
      INTEGER ISEVOK(*)
      INTEGER NUMLNK(*), NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*(NAMLEN) ATNAME(*)
      INTEGER DBMAPEL(*)
      INTEGER DBMAPND(*)
      INTEGER MAPND(*)
      INTEGER MAPEL(*)
      INTEGER IDNPS(*), NNNPS(*), IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      INTEGER IDESS(*), NEESS(*), NNESS(*), IXEESS(*), IXNESS(*)
      INTEGER LTEESS(*), LTSESS(*)
      REAL FACESS(*)
      REAL VARGL(*), VARNP(*), VAREL(*)

      integer N, L, KICHECK, KISCR, KRCHECK, NERR, MEM, NSTEP
      REAL TIME

      L = MAX (NUMEL, NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL, NUMNP)
      CALL MDRSRV ('ICHECK', KICHECK, L)
      CALL MDRSRV ('ISCR',   KISCR, LESSEL)
      CALL MDRSRV ('RCHECK', KRCHECK, NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 240

      write (*,*)
      CALL CKMAP (NUMNP, DBMAPND, IA(KICHECK), 'Node   ')
      CALL CKMAP (NUMEL, DBMAPEL, IA(KICHECK), 'Element')
      CALL CKELB (NELBLK, NUMEL, NUMNP, EBTYPE, 
     &  IDELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB, ATNAME, 
     &  IA(KICHECK), MAPND)
      CALL CKNPS (NUMNPS, LNPSNL, NUMNP,
     &  IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, A(KICHECK))
      CALL CKESS (NUMESS, LESSEL, LESSNL, NUMEL, NUMNP,
     &  IDESS, NEESS, NNESS, IXEESS, IXNESS,
     &  LTEESS, LTSESS, FACESS,
     *  A(KISCR), A(KICHECK), A(KRCHECK), NDIM,
     *  MAPEL, MAPND)
      
 240  CONTINUE
      CALL MDDEL ('ICHECK')
      CALL MDDEL ('RCHECK')
      
      IF (EXODUS) THEN
        DO 250 N = 1, NSTEPS
          NSTEP = N
          CALL TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK,
     &      TIME, VARGL, VARNP, VAREL)
          IF (N .NE. NSTEP) GOTO 260
 250    CONTINUE
 260    CONTINUE
        
        NSTEP = 1
        CALL TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK,
     &    TIME, VARGL, VARNP, VAREL)
      END IF
      
      WRITE (*, *)
      WRITE (*, 10000) 'Database check is completed'
      
10000 FORMAT (1X, 5A)
      return
      end
      
