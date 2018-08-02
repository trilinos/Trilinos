C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Id: loctol.f,v 1.1 1991/02/21 15:44:03 gdsjaar Exp $
C $Log: loctol.f,v $
C Revision 1.1  1991/02/21 15:44:03  gdsjaar
C Initial revision
C
      SUBROUTINE LOCTOL (TYPE, NDIM, RV, KV)
C
C     This routine is used to set the tolerances and distances
C        used in the LOCATE routines.
C     If a tolerance is not entered (blank field), then
C        the tolerance is set to the entered distance value, and
C        the distance is set to 0.0
C     If a tolerance is entered, the values are returned with no
C        changes
C
      DIMENSION RV(*), KV(*)
      CHARACTER*(*) TYPE
      LOGICAL MATSTR

         IF (MATSTR(TYPE, 'LINE', 1)) THEN
            IF (NDIM .EQ. 3) THEN
               IF (KV(11) .EQ. -1) THEN
                  RV(11) = RV(10)
                  RV(10) = 0.0
               END IF
            ELSE
               IF (KV(9) .EQ. -1) THEN
                  RV(9) = RV(8)
                  RV(8) = 0.0
               END IF
            END IF
         ELSE IF (MATSTR(TYPE, 'PLANE', 2)) THEN
            IF (NDIM .EQ. 3) THEN
               IF (KV(11) .EQ. -1) THEN
                  RV(11) = RV(10)
                  RV(10) = 0.0
               END IF
            ELSE
               CONTINUE
            END IF
         ELSE IF (MATSTR(TYPE, 'POINT', 2)) THEN
            IF (NDIM .EQ. 3) THEN
               IF (KV(8) .EQ. -1) THEN
                  RV(8) = RV(7)
                  RV(7) = 0.0
               END IF
            ELSE
               IF (KV(7) .EQ. -1) THEN
                  RV(7) = RV(6)
                  RV(6) = 0.0
               END IF
            END IF
         END IF
      RETURN
      END
