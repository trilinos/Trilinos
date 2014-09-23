C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: dread.f,v 1.2 1998/07/14 18:18:42 gdsjaar Exp $
C $Log: dread.f,v $
C Revision 1.2  1998/07/14 18:18:42  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:06:24  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:06:21  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]DREAD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)
C***********************************************************************
C
C  SUBROUTINE DREAD = SETS ALL PARAMETERS UP FOR READING FROM A DIGI-PAD
C
C***********************************************************************
C
C  SUBROUTINES CALLED:
C     DPREAD   = READS INPUT FROM A DIGI-PAD DIGITIZER
C
C***********************************************************************
C
C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C     BUTTON = THE MOUSE BUTTON PUSHED
C     X1     = THE X ZERO SHIFT IN USER COORDINATES
C     Y1     = THE Y ZERO SHIFT IN USER COORDINATES
C     XX1    = THE X ZERO SHIFT IN DIGITIZED COORDINATES
C     YY1    = THE Y ZERO SHIFT IN DIGITIZED COORDINATES
C     SCALE  = THE SCALE FACTOR FROM DIGITIZED TO USER COORDINATES
C     CT     = THE COSINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     ST     = THE SINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     XNEW   = THE NEW DIGITIZED X VALUE BEFORE TRANSFORMATIONS
C     YNEW   = THE NEW DIGITIZED Y VALUE BEFORE TRANSFORMATIONS
C     NCB    = THE NUMBER OF BUTTONS ON THE MOUSE  (BIT-PAD-ONE)
C     DEL    = THE DELTA DISTANCE BETWEEN ACCEPTABLE POINTS  (TALOS)
C
C***********************************************************************
C
      CHARACTER*1 BUTTON
C
      CALL DPREAD (XNEW, YNEW, BUTTON)
      XNEW = XNEW - XX1
      YNEW = YNEW - YY1
      X = X1 + SCALE * ( (CT * XNEW) + (ST * YNEW))
      Y = Y1 + SCALE * ( (-ST * XNEW) + (CT * YNEW))
      RETURN
      END
