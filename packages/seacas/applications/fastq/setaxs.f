C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: setaxs.f,v 1.1 1990/11/30 11:15:21 gdsjaar Exp $
C $Log: setaxs.f,v $
C Revision 1.1  1990/11/30 11:15:21  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]SETAXS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SETAXS (X, Y)
C***********************************************************************
C
C  SETAXS = SETS UP THE AXIS AS NEEDED FOR PLOTTING
C
C***********************************************************************
C
      DIMENSION BUFF (11), X (2), Y (2)
C
C  GET THE AXIS ON THE CORRECT LOCATION ON THE GRAPH
C
      CALL MP2PT (1, X (1), Y (1), X01, Y01, MASK)
      CALL MP2PT (1, X (2), Y (2), X02, Y02, MASK)
C
C  FORCE X AND Y LIMITS ON THE GRAPH
C
      CALL PLTSTG (1, X01)
      CALL PLTSTG (2, Y01)
      CALL PLTSTG (3, X02 - X01)
      CALL PLTSTG (4, Y02 - Y01)
C
C  TURN OFF THE ZERO LINE PLOT
C
      CALL PLTSTG (37, 0.)
C
C  GET NICE INTERVALS ON THE AXIS
C
      CALL PLTINI (X (1), X (2), XSTART, XEND, XINT, IXEXP, IXTIC)
      CALL PLTINI (Y (1), Y (2), YSTART, YEND, YINT, IYEXP, IYTIC)
C
C  SET ALL THE BUFFER PARAMETERS
C
      BUFF (1) = 4.
      BUFF (2) = X (1)
      IF (IXEXP .EQ. 0) THEN
         BUFF (3) = XSTART
         BUFF (5) = XINT
      ELSE
         BUFF (3) = XSTART *  (10. ** FLOAT (IXEXP))
         BUFF (5) = XINT *  (10. ** FLOAT (IXEXP))
      ENDIF
      BUFF (4) = X (2)
      BUFF (6) = 1.
      BUFF (7) = Y (1)
      IF (IYEXP .EQ. 0) THEN
         BUFF (8) = YSTART
         BUFF (10) = YINT
      ELSE
         BUFF (8) = YSTART *  (10. ** FLOAT (IYEXP))
         BUFF (10) = YINT *  (10. ** FLOAT (IYEXP))
      ENDIF
      BUFF (9) = Y (2)
      BUFF (11) = 1
C
C  FORCE THE CORRECT AXIS SETUP
C
      CALL PLTSTG (11, BUFF)
C
C  PLOT THE AXIS
C
      CALL PLTGPH (X, Y,  - 2, 'X', ' ', 'Y', ' ')
C
C  PUT THE CLIPPING RECTANGLE RIGHT AT THE AXIS
C
      CALL MP2PT (1, X (1), Y (1), X01, Y01, MASK)
      CALL MP2PT (1, X (2), Y (2), X02, Y02, MASK)
      CALL MPVIEW (X01, X02, Y01, Y02)
      CALL MPORT2 (X (1), X (2), Y (1), Y (2))
C
      RETURN
C
      END
