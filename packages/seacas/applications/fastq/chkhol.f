C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
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

C $Id: chkhol.f,v 1.3 1999/06/21 22:43:40 gdsjaar Exp $
C $Log: chkhol.f,v $
C Revision 1.3  1999/06/21 22:43:40  gdsjaar
C Fixed more uninitialized variables; one was causing core dump on g77
C compiled executable.
C
C VERSN was not consistently defined -- now 10 characters everywhere
C
C Updated so full version string output
C
C Added capability to debug memory using unit specified in EXT99
C variable. Similar to STRTUP in SUPLIB
C
C Cleaned up some other code
C
C Upped version
C
C Revision 1.2  1998/07/14 18:18:27  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:04:31  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:04:30  gdsjaar
c Initial revision
c 
CC* FILE: [.QMESH]CHKHOL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO CHKHOL TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE CHKHOL (IA, L, MP, ML, MS, MR, MSC, IPOINT, COOR,
     &   IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN,
     &   ISIDE, NLPS, IFLINE, ILLIST, IREGN, NSPR, IFSIDE, ISLIST,
     &   NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB,
     &   IFHOLE, NHPR, IHLIST, LINKP, LINKL, LINKS, LINKR, LINKSC,
     &   LINKPB, LINKLB, LINKSB, RSIZE, NPREGN, NPSBC, NPNODE, MAXNP,
     &   MAXNL, MXNPER, MXRNBC, MXRSBC, X, Y, NID, LISTL, MARKED, MXNL,
     &   MAXNBC, MAXSBC, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG,
     &   LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH,
     &   REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &   NOROOM, ERRCHK, ERR)
C***********************************************************************
C
C  CHKRGN - CHECK THAT A REGION MAY BE MESHED
C
C***********************************************************************
C
      DIMENSION IA(1)
      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML), ISIDE(MS), NLPS(MS)
      DIMENSION IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION RSIZE (MR)
      DIMENSION NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKR(2, MR), LINKSC(2, MR), LINKPB(2, MP)
      DIMENSION LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION X(MAXNP), Y(MAXNP), NID(MAXNP)
      DIMENSION LISTL(MAXNL), MARKED(3, MAXNL)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
C
      DIMENSION IDUMMY(1)
C
      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL NOROOM, EVEN, ERR, CCW, REAL, ADDLNK, REMESH
      LOGICAL COUNT, ERRCHK
C
      addlnk = .false.
      COUNT = .TRUE.
      EVEN = .FALSE.
      REAL = .FALSE.
C
C  CHECK TO MAKE SURE CONNECTING DATA FOR THE REGION EXISTS
C  AND FILL IN ANY BLANK INTERVALS ACCORDING TO THE GIVEN SIZE
C  FOR THE REGION AND THE LINE'S LENGTH
C
      IF (NHPR(L) .GT. 0) THEN
         DO 100 I = IFHOLE(L), IFHOLE(L) + NHPR(L) - 1
            IPNTR1 = 0
            CALL LTSORT (MR, LINKR, IHLIST(I), IPNTR1, ADDLNK)
            IF (IPNTR1 .GT. 0) THEN
               LL = IPNTR1
               CALL DATAOK (MP, ML, MS, MR, LL, IREGN(LL), COOR, ILINE,
     &            LTYPE, NINT, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, RSIZE(LL), ERRCHK, ERR)
               IF (ERR) THEN
                  WRITE (*, 10000) IREGN(LL)
                  ADDLNK = .TRUE.
                  IMINUS = -LL
                  CALL LTSORT (MR, LINKR, IREGN(LL), IMINUS, ADDLNK)
                  ADDLNK = .FALSE.
C
C  CALCULATE THE PERIMETER OF THE REGION
C
               ELSE
                  KNBC = 0
                  KSBC = 0
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO PERIM TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
                  CALL PERIM (MP, ML, MS, NSPR(LL), MAXNL, MAXNP, 1, 1,
     &               KNBC, KSBC, IREGN(LL), IPOINT, COOR, IPBOUN, ILINE,
     &               LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE,
     &               NLPS, IFLINE, ILLIST, ISLIST(IFSIDE(LL)), NPPF,
     &               IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB,
     &               LISTSB, LINKP, LINKL, LINKS, LINKPB, LINKLB,
     &               LINKSB, X, Y, NID, NPER, LISTL, NL, IDUMMY,
     &               MARKED, EVEN,  REAL, ERR, CCW, COUNT, NOROOM,
     &               AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG,
     &               LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD,
     &               NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &               IDIVIS, SIZMIN, EMAX, EMIN)
                  IF ((NPER .LE. 0) .OR. (ERR)) THEN
                     WRITE (*, 10010) IREGN(LL)
                     ADDLNK = .TRUE.
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/23/90
CC* MODIFICATION: PUT THE CORRECT POINTER INTO THE HOLE REGION LINK SLOT
C
                     IMINUS = -LL
                     CALL LTSORT (MR, LINKR, IREGN(LL), IMINUS, ADDLNK)
                     ADDLNK = .FALSE.
                  ELSE
C
C  WHEN CHECKING THE MAXIMUMS - ADD ENOUGH FOR ONE MORE INTERVAL
C  ON THE LINE AS THIS LINE MAY BE INCREMENTED BY ONE IF THE
C  PERIMETER IS ODD
C
                     MAXNBC = MAX(MAXNBC, KNBC + 3 + MXRNBC)
                     MAXSBC = MAX(MAXSBC, KSBC + 3 + MXRSBC)
                     MXNL   = MAX(MXNL, NL)
                     MXNPER = MAX(MXNPER, NPER + 2)
C
C  MARK THE LINES AND POINTS IN THE REGION AS BEING USED
C
                     CALL MKUSED (MAXNL, MP, ML, LISTL, IPOINT, NINT,
     &                  LINKP, LINKL, LCON, NL)
                  ENDIF
               ENDIF
            ELSE
               WRITE (*, 10020) IREGN(LL)
               ERR = .TRUE.
            ENDIF
  100    CONTINUE
      END IF
C
      RETURN
C
10000 FORMAT (' ** ERROR - DATA PROBLEMS FOR HOLE REGION:', I5, ' **')
10010 FORMAT (' ** ERROR - PERIMETER GENERATION ERRORS FOR HOLE REGION:'
     &   , I5, ' **')
10020 FORMAT (' ** ERROR - HOLE REGION', I5, ' DOES NOT EXIST **')
      END
