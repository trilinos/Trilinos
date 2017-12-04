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

C $Id: add2cn.f,v 1.1 1990/11/30 11:02:40 gdsjaar Exp $
C $Log: add2cn.f,v $
C Revision 1.1  1990/11/30 11:02:40  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADD2CN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADD2CN TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE ADD2CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &   ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, NLOOP,
     &   I1, IAVAIL, NAVAIL, GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD,
     &   YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD,
     &   NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN,
     &   EMAX, EMIN)
C***********************************************************************
C
C  SUBROUTINE ADD2CN = ADDS TWO CENTER NODES IN THE MIDDLE OF 6 CLOSING
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL GRAPH, VIDEO, ERR, AMBIG, SIZEIT, NOROOM
C
      ERR = .FALSE.
      AMBIG = .FALSE.
C
      I2 = LNODES (3, I1)
      I3 = LNODES (3, I2)
      I4 = LNODES (3, I3)
      I5 = LNODES (3, I4)
      I6 = LNODES (3, I5)
C
C  CALCULATE THE TWO NEW CENTER LOCATIONS
C
      X16 = ( XN(I1) + XN(I6) ) * .5
      X45 = ( XN(I4) + XN(I5) ) * .5
      XMID = (X16 + X45) * .5
      XNEW1 = (X16 + XMID) * .5
      XNEW2 = (X45 + XMID) * .5
C
      Y16 = ( YN(I1) + YN(I6) ) * .5
      Y45 = ( YN(I4) + YN(I5) ) * .5
      YMID = (Y16 + Y45) * .5
      YNEW1 = (Y16 + YMID) * .5
      YNEW2 = (Y45 + YMID) * .5
C
      ZERO = 0.
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADDNOD TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      CALL ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN, ANGLE, BNSIZE,
     &   LNODES, XNEW1, YNEW1, ZERO, NNN, KKK, LLL, I6, I1, I2,
     &   AMBIG, IDUM, SIZEIT, ERR, NOROOM, XNOLD, YNOLD, NXKOLD,
     &   LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &   REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
      IF ((ERR) .OR. (NOROOM)) GOTO 100
      INEW = NNN
      LNODES (4, INEW) = - 2
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL D2NODE (MXND, XN, YN, I6, INEW)
         CALL D2NODE (MXND, XN, YN, I2, INEW)
         CALL SFLUSH
         IF (VIDEO) CALL SNAPIT (1)
      ENDIF
      CALL ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN, ANGLE, BNSIZE,
     &   LNODES, XNEW2, YNEW2, ZERO, NNN, KKK, LLL, INEW, I2, I3,
     &   AMBIG, IDUM, SIZEIT, ERR, NOROOM, XNOLD, YNOLD, NXKOLD,
     &   LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &   REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
      IF ((ERR) .OR. (NOROOM)) GOTO 100
      JNEW = NNN
      LNODES (4, JNEW) = - 2
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL D2NODE (MXND, XN, YN, INEW, JNEW)
         CALL D2NODE (MXND, XN, YN, I3, JNEW)
         CALL SFLUSH
         IF (VIDEO) CALL SNAPIT (1)
      ENDIF
      CALL CONNOD (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN, ANGLE,
     &   LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, JNEW, I3, I4, I5, IDUM,
     &   NLOOP, IAVAIL, NAVAIL, GRAPH, VIDEO, NOROOM, ERR)
      IF ((NOROOM) .OR. (ERR)) GOTO 100
      CALL CLOSE4 (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   INEW, JNEW, I5, I6, KKK, ERR)
C
  100 CONTINUE
C
      RETURN
C
      END
