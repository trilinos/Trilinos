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

C $Log: upcase.f,v $
C Revision 1.2  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:17:19  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:59:25  gdsjaar
c Added RCS Id and Log to all files
c
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c
c ROUTINE:              Upcase
c
c DESCRIPTION:          Changes ASCII text strings to upper case.
c
c AUTHOR:               John H. Glick
c                       Sandia National Laboratories
c                       Division 1511
c
c DATE:                 December 20, 1988
c
c TYPE OF SUBPROGRAM:   subroutine
c
c USEAGE:               call upcase ( string )
c
c PARAMETERS:
c
c        character*(*) string -- ( Input and Output )
c                       String that is to be converted to
c                          upper case.
c
c CALLS:
c
c        len (INTRINSIC) -- returns length of character string.
c        ichar (INTRINSIC) -- returns ASCII integer value
c                             of passed character.
c        char (INTRINSIC) -- returns ASCII character assigned
c                            to passed integer value.
c
c
c GLOBAL VARIABLES REFERENCED:   none
c
c SYSTEM DEPENDENCIES:           none
c
c CALLING ROUTINE(S):            filhnd (BLOT)
c
c ======================================================================
c ======================================================================
c
      subroutine upcase ( string )
c
c **********************************************************************
c
c        parameter
c
      character*(*) string
c
c **********************************************************************
c
c        declarations
c
      integer length       
c          length of character string
      integer ccode        
c          integer id of an ASCII character
c
c***********************************************************************
c***********************************************************************
c
      length = len(string)
      do 100 i = 1,length
         ccode = ichar(string(i:i))
         if ( ccode .ge. 97  .AND.  ccode .le. 122 )
     &      string(i:i) = char(ccode-32)
  100 continue

      return
      end
