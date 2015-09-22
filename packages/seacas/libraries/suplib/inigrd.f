C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C=======================================================================
      SUBROUTINE INIGRD (FROM, TO, GRAD, NINT, NPTS, ARRAY)
C=======================================================================

C   $Id: inigrd.f,v 1.4 2009/03/25 12:46:02 gdsjaar Exp $
C   $Log: inigrd.f,v $
C   Revision 1.4  2009/03/25 12:46:02  gdsjaar
C   Add copyright and license notice to all files.
C
C   Revision 1.3  1992/01/02 18:55:16  gdsjaar
C   Moved from gen3d to here, overwrite old one
C
c Revision 1.2  1991/02/07  10:35:19  gdsjaar
c Removed Debug Comment lines from files
c
c Revision 1.1.1.1  90/08/20  12:21:56  gdsjaar
c Gen3D Mesh Generation Program
c 
c Revision 1.1  90/08/20  12:21:55  gdsjaar
c Initial revision
c 
C     INIGRD: Initialize array ARRAY with NPTS values.
C             values are calculated to have a gradient of GRAD
C             ranging from FROM to TO with NINT segments.
C
C --- FROM - IN - Minimum value of range
C --- TO   - IN - Maximum value of range
C --- GRAD - IN - Gradient
C --- NINT - IN - Number of intervals in range
C --- NPTS - IN - Number of points in range.  Can be less than NINT if
C                 do not need full range, but want gradient spacing
C                 based on full range.  Normally NPTS = NINT + 1
C --- ARRAY- OUT- Range of values
C
      REAL    ARRAY(NPTS)
      LOGICAL NOGRAD
C
      IF (FROM .EQ. TO  .OR. NINT .LE. 0) THEN
         CALL PRTERR ('PROGRAM',
     *      'invalid values passed to INIGRD')
         RETURN
      END IF

      NOGRAD = (ABS(GRAD - 1.0)/NINT .LE. 1.0e-7)
C
      IF (NOGRAD) THEN
         D3 = 1.0 / NINT
         DO 10 I=1, NPTS
            ARRAY(I) = D3 * (I-1)
   10    CONTINUE
      ELSE
         ARRAY(1) = 0.0
         D3 = (1.0 - GRAD) / (1.0 - GRAD**NINT)
         DO 20 I=2,NPTS
            ARRAY(I) = ARRAY(I-1) + D3 * GRAD**(I-2)
   20    CONTINUE
      END IF

      DO 30 I=1, NPTS
         ARRAY(I) = FROM + ARRAY(I) * (TO - FROM)
   30 CONTINUE
C
      RETURN
      END
