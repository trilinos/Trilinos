C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C * Neither the name of Sandia Corporation nor the names of its
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
C 

C   $Id: mirss.f,v 1.3 1999/01/27 15:08:25 gdsjaar Exp $
C=======================================================================
      SUBROUTINE MIRSS (IDFRO, IDBCK,
     &   NESUR, NESFRO, NESBCK, LTSSS3)
C=======================================================================

C   --*** MIRSS *** (GEN3D) Modifies sideset node order to account
C                           for mirroring about axes
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from WRESS written by Amy Gilkey
C   --
C   --Parameters:
C   --   IDFRO  - IN - ids for front surface side sets; (0) = length
C   --   IDBCK  - IN - ids for back surface side sets; (0) = length
C   --   NESUR  - IN - the number of elements in the surface side set
C   --   NESFRO - IN - the elements in the front surface side set
C   --   NESBCK - IN - the elements in the back surface side set
C   --   LTEES3 - IN - the element faces for all 3D sets
C   --
C   --Common Variables:
C   --   Uses NUMESS of /DBNUMS/
C   --   Uses LESSEO of /DBNUM3/
C   --
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER NESFRO(*), NESBCK(*)
      INTEGER LTSSS3(*)

      LOGICAL ANYESS

      INTEGER NEWFAC(6)
      DATA NEWFAC /4, 3, 2, 1, 5, 6/

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYESS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMESS .GT. 0)

C   --Write 3D

      IF (ANYESS) THEN
         DO 10 NL = 1, LESSEO
C ... non-front and non-back sidesets
           LTSSS3(NL) = NEWFAC(LTSSS3(NL))
   10    CONTINUE

C ... Front and back don't get mirrored...(?)
         DO 20 NL = 1, NESUR
   20    CONTINUE
         DO 30 NL = 1, NESUR
   30    CONTINUE
      END IF

      RETURN
      END
