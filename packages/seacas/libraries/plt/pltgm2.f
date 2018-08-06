C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Id: pltgm2.f,v 1.1 1993/07/16 16:48:14 gdsjaar Exp $ 
C $Log: pltgm2.f,v $
C Revision 1.1  1993/07/16 16:48:14  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTGM2(XL,XU,YL,YU,PXL,PXU,PYL,PYU,UMAP)
      REAL UMAP(*)

      DU = XU - XL
      DP = PXU - PXL
      IF (DU.EQ.0.) THEN
         DU = 1.
      END IF

      UMAP(1) = DP/DU
      UMAP(1+1) = 0.
      UMAP(1+2) = 0.
      DV = YU - YL
      DQ = PYU - PYL
      IF (DV.EQ.0.) THEN
         DV = 1.
      END IF

      UMAP(1+3) = DQ/DV
      UMAP(5) = PXL - XL*UMAP(1)
      UMAP(5+1) = PYL - YL*UMAP(1+3)
      UMAP(7) = PXL
      UMAP(7+1) = PYL
      UMAP(9) = PXU
      UMAP(9+1) = PYL
      UMAP(11) = PXU
      UMAP(11+1) = PYU
      UMAP(13) = PXL
      UMAP(13+1) = PYU
      RETURN

      END
