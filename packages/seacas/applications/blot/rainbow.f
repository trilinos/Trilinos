      subroutine rainbow(h, s, v, r, g, b) 
      real h, s, v
      real r, g, b
      
C This routine computes colors suitable for use in color level plots.
C Typically s=v=1 and h varies from 0 (red) to 1 (blue) in
C equally spaced steps.  (h=.5 gives green; 1<h<1.5 gives magenta.)
C To convert for frame buffer, use   R = floor(255.999*pow(*r,1/gamma))  etc.
C To get tables calibrated for other devices or to report complaints,
C contact ehg@research.att.com.

C  The author of this software is Eric Grosse.  Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C  The author of this software is Eric Grosse.  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C  The author of this software is Eric Grosse.  certain rights in this software.
C  The author of this software is Eric Grosse.          
C  The author of this software is Eric Grosse.  Redistribution and use in source and binary forms, with or without
C  The author of this software is Eric Grosse.  modification, are permitted provided that the following conditions are
C  The author of this software is Eric Grosse.  met:
C  The author of this software is Eric Grosse.  
C  The author of this software is Eric Grosse.      * Redistributions of source code must retain the above copyright
C  The author of this software is Eric Grosse.        notice, this list of conditions and the following disclaimer.
C  The author of this software is Eric Grosse.  
C  The author of this software is Eric Grosse.      * Redistributions in binary form must reproduce the above
C  The author of this software is Eric Grosse.        copyright notice, this list of conditions and the following
C  The author of this software is Eric Grosse.        disclaimer in the documentation and/or other materials provided
C  The author of this software is Eric Grosse.        with the distribution.
C  The author of this software is Eric Grosse.      * Neither the name of Sandia Corporation nor the names of its
C  The author of this software is Eric Grosse.        contributors may be used to endorse or promote products derived
C  The author of this software is Eric Grosse.        from this software without specific prior written permission.
C  The author of this software is Eric Grosse.  
C  The author of this software is Eric Grosse.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C  The author of this software is Eric Grosse.  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C  The author of this software is Eric Grosse.  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C  The author of this software is Eric Grosse.  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C  The author of this software is Eric Grosse.  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C  The author of this software is Eric Grosse.  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C  The author of this software is Eric Grosse.  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C  The author of this software is Eric Grosse.  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C  The author of this software is Eric Grosse.  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C  The author of this software is Eric Grosse.  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C  The author of this software is Eric Grosse.  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C  Permission to use, copy, modify, and distribute this software for any
C  purpose without fee is hereby granted, provided that this entire notice
C  is included in all copies of any software which is or includes a copy
C  or modification of this software and in all copies of the supporting
C  documentation for such software.
C  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
C  WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
C  REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
C  OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.

      INTEGER i
      REAL huettab(0:60)
      DATA huettab /0.0000, 0.0062, 0.0130, 0.0202, 0.0280, 
     *              0.0365, 0.0457, 0.0559, 0.0671, 0.0796,
     *              0.0936, 0.1095, 0.1275, 0.1482, 0.1806,
     *              0.2113, 0.2393, 0.2652, 0.2892, 0.3119,
     *              0.3333, 0.3556, 0.3815, 0.4129, 0.4526,
     *              0.5060, 0.5296, 0.5501, 0.5679, 0.5834,
     *              0.5970, 0.6088, 0.6191, 0.6281, 0.6361,
     *              0.6430, 0.6490, 0.6544, 0.6590, 0.6631,
     *              0.6667, 0.6713, 0.6763, 0.6815, 0.6873,
     *              0.6937, 0.7009, 0.7092, 0.7190, 0.7308,
     *              0.7452, 0.7631, 0.7856, 0.8142, 0.8621,
     *              0.9029, 0.9344, 0.9580, 0.9755, 0.9889, 1.0000/

C computed from the FMC-1 color difference formula
C Barco monitor, max(r,g,b)=1, n=61 magenta,  2 Jan 1986

      
      H = 60.0 * MOD(H / 1.5, 1.)
      I = INT(H)
      H = huettab(i) + (huettab(i+1) - huettab(i)) * (h - i)
      CALL dhsv2rgb(h, s, v, r, g, b)
      RETURN
      END

      SUBROUTINE dhsv2rgb(h, s, v, r, g, b)
C...hexcone model...
      REAL h, s, v
      REAL r, g, b
C...all variables in range [0,1[ 
C...here, h=.667 gives blue, h=0 or 1 gives red.
C...see Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78

      INTEGER i
      REAL f, m, n, k

      h = 6 * mod(h, 1.)
      i = int(h)
      f = h - i
      m = (1 - s)
      n = (1 - s * f)
      k = (1 - (s * (1 - f)))

      if (i .eq. 0) then
         r = 1 
         g = k 
         b = m 
      else if (i .eq. 1) then
         r = n 
         g = 1 
         b = m 
      else if (i .eq. 2) then
         r = m 
         g = 1 
         b = k 
      else if (i .eq. 3) then
         r = m 
         g = n 
         b = 1 
      else if (i .eq. 4) then
         r = k 
         g = m 
         b = 1 
      else if (i .eq. 5) then
         r = 1 
         g = m 
         b = n 
      end if
      f = max(r, g, b)
      f = v / f
      r = r * f
      g = g * f
      b = b * f
        
      return
      end
