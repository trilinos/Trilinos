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
      BLOCK DATA BLKDAT
C=======================================================================

C   --*** BLKDAT *** (BLOT) Block data
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --BLKDAT contains all the COMMON areas defined for BLOT.

C=======================================================================
C                               B L O T
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
C      --These parameters define the indices of 2D and 3D limit arrays

      include 'params.blk'
      include 'progqa.blk'
      include 'dbase.blk'
      include 'dbtitl.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'dbnams.blk'
      include 'legopt.blk'
      include 'times.blk'
      include 'layout.blk'
      include 'plcolr.blk'
      include 'cmap-lst.blk'
      include 'light.blk'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

C=======================================================================
C                        M E S H   P L O T

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)
C      --These parameters define the mesh display (see MSHLIN of /MSHOPT/)

      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'cutopt.blk'
      include 'layoud.blk'
      include 'devdat.blk'
      include 'pick.blk'
      include 'sizes.blk'
      include 'linthc.blk'
      include 'sphele.blk'
C
C=======================================================================
C                             D E T O U R

      include 'detopt.blk'
      include 'cntr.blk'
      include 'icexct.blk'
      include 'icrnbw.blk'
      include 'etcopt.blk'
      include 'nodzom.blk'
      include 'axsplt.blk'
C=======================================================================
C                          X - Y   P L O T

      include 'xyopt.blk'
      include 'xylim.blk'
      include 'xylab.blk'
      include 'neutr.blk'

C=======================================================================
C                        P A T H L I N E

      include 'lnvars.blk'
C=======================================================================
C                             T P L O T


      include 'tpvars.blk'
C=======================================================================
C                             S P L O T

      include 'selne.blk'
      include 'spvars.blk'

C=======================================================================

      DATA IEXCON /0/
      END
