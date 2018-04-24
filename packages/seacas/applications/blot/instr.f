C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Log: instr.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:53  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:43  gdsjaar
c Added RCS Id and Log to all files
c
c ======================================================================
c ======================================================================
c ======================================================================
c ======================================================================
c
c ROUTINE:              instr
c
c DESCRIPTION:          Returns .TRUE. if a command file was specified
c                       on the BLOT command line, indicating
c                       that instructions are to be read from that file.
c                       Returns .FALSE. if no file was specified on the
c                       BLOT command line, indicating
c                       that instructions are to be read interactively.
c
c AUTHOR:               John H. Glick
c                       Sandia-2017 National Laboratories
c                       Division 1511
c
c DATE:                 December 20, 1988
c
c TYPE OF SUBPROGRAM:   logical function
c
c USEAGE:               instr()
c
c PARAMETERS:           none
c
c CALLS:                exname ( SUPES )
c
c GLOBAL VARIABLES REFERENCED:   none
c
c CALLING ROUTINE(S):   getins (BLOT)
c
c SYSTEM DEPENDENCIES:  none
c
c ======================================================================
c ======================================================================
c
      logical function instr()

c
c **********************************************************************
c
c        declarations
c
      integer symnum       
c          the number of the system command symbol defining if a 
c          command file was specified.
      character*12 name    
c          string that holds the value of the
c          symbol associated with the symnum
c          = "YES" if a file was specified
c            in the P4 field of the TRINITY
c            command line, indicating that
c            instructions are to be read
c            from that file
c          = "NO" otherwise.
      integer ln           
c          length of the string stored in name
c
c **********************************************************************
c
c        data statements
c
      data symnum / 4 /
c **********************************************************************
c **********************************************************************
c

c
c           get symbol from operating system
c
      call exname ( -symnum, name, ln )

c
c           print symbol for debugging purposes
c
c     print *, 'name = ',name,' ln = ',ln
c
      if ( name(1:ln) .eq. 'YES' ) then
         instr = .TRUE.
      else
         instr = .FALSE.
      endif

      return
      end
