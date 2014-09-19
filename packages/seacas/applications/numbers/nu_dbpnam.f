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

C=======================================================================
      SUBROUTINE DBPNAM (OPTION, NVARGL, NVARNP, NVAREL,
     &                   NAMEGV, NAMENV, NAMEEV)
C=======================================================================
C$Id: dbpnam.f,v 1.1 1999/02/16 21:37:59 gdsjaar Exp $

C   --*** DBPNAM *** (EXOLIB) Print database variable names
C   --   Written by Amy Gilkey - revised 01/21/88
C   --
C   --DBPNAM displays the database variable names.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --                 'G' to print global variable names
C   --                 'N' to print nodal variable names
C   --                 'E' to print element variable names
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)
C   --   NAMEGV - IN - the global variable names (if OPTION)
C   --   NAMENV - IN - the nodal variable names (if OPTION)
C   --   NAMEEV - IN - the element variable names (if OPTION)

      include 'exodusII.inc'

      CHARACTER*(*) OPTION
      INTEGER NVARGL, NVARNP, NVAREL
      CHARACTER*(MXSTLN) NAMEGV(*)
      CHARACTER*(MXSTLN) NAMENV(*)
      CHARACTER*(MXSTLN) NAMEEV(*)
      LOGICAL ALL

      ALL = (OPTION .EQ. '*')

      WRITE (*, 10000)

      IF (ALL .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
         WRITE (*, 10010) 'Global: ', (NAMEGV(I), I=1,NVARGL)
      END IF
      IF (ALL .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         WRITE (*, 10010) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
      END IF
      IF (ALL .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         WRITE (*, 10010) 'Element:', (NAMEEV(I), I=1,NVAREL)
      END IF

      RETURN

10000  FORMAT (/, 1X, 'Variables Names:')

10010  FORMAT (4X, A, :, 2 (2X, A32), :, /,
     &        (12X, 2 (2X, A32)))
      END
