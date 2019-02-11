C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
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

C=======================================================================
      SUBROUTINE LNKFNC (NUMSTO, *)
C=======================================================================

C   --*** LNKFNC *** (ALGEBRA) Assign storage for time functions
C   --   Written by Amy Gilkey - revised 07/22/87
C   --
C   --LNKFNC sets up the storage locations for the time functions that
C   --need storage for results that must be saved over time steps.
C   --
C   --Parameters:
C   --   NUMSTO - IN/OUT - the number of variable storage locations needed
C   --   * - return statement if an error is found; message is printed
C   --
C   --Common Variables:
C   --   Sets ITMENT of /ENT../
C   --   Uses NUMEQN, NUMENT, TYPENT, INXENT of /ENT../
C   --   Uses FNCSTO of /FNCTB./

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_ent.blk'
      include 'ag_fnctbc.blk'

C   --Allocate storage for time functions

      DO 110 NEQN = 1, NUMEQN
         DO 100 NENT = 3, NUMENT(NEQN)
            IF (TYPENT(NENT,NEQN) .EQ. 'F') THEN
               INX = INXENT(NENT,NEQN)
               IF (FNCSTO(INX)) THEN
                  NUMSTO = NUMSTO + 1
                  ITMENT(NENT,NEQN) = NUMSTO
               ELSE
                  ITMENT(NENT,NEQN) = 0
               END IF
            END IF
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
