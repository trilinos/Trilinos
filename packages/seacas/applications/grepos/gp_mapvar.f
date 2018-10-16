C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C * Neither the name of NTESS nor the names of its
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

C=======================================================================
      subroutine mapvar(nold, nnew, nvar, map, vars, scr)
C=======================================================================
C   --*** MAPVAR *** (GREPOS) Map values from old array to new.
C   --
C   --MAPVAR reorders the data in VARS based on the map in MAP
C   --
C   --Parameters:
C   --   NOLD - IN - number of old values, 
C   --   NNEW - IN - number of new values
C   --   NVAR - IN - the number of variables
C   --   NVARDM - IN - the row dimension of VARNP
C   --   MAP  - IN - map from new value to old MAP(NEW) = OLD
C                    size is 'nnew'      
C   --   VARS - IN/OUT - the values. On input in old order,
C                        on output in new order      
C   --   SCR  - TMP - temporary storage area

      
      integer map(*)
      real vars(*)
c     real vars(nold, nvar)
      real scr(*)

C ... VARS should really be a doubly-dimensioned array (NOLD, NVAR),
C     The dimensions need to be in this order so we can read them
C     in using exgev and exgnv.  But, this order doesn't work very
C     well when trying to reduce the size of NOLD

C ... TODO: Need to use the truth table to make sure variables
C           exist for each element.      
      do 30 ivar = 1, nvar
        do 10 i = 1, nnew
          scr(i) = vars(map(i) + nold * (ivar-1) )
 10     continue
        
        do 20 i = 1, nnew
          vars(i + nnew * (ivar-1)) = scr(i)
 20     continue
 30   continue

      return
      end
