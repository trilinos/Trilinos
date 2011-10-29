C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C 

C -*- Mode: fortran -*-
C=======================================================================
C $Id: offset.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: offset.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.3  1997/11/18 15:55:30  gdsjaar
C Changes to improve efficiency (factor of 2 to 15 on some Goodyear Treads)
C
C Redid nodeset munching to reduce number of locint calls. First do a
C quick scan on array to get maximum node id in searched array. Then,
C before doing 'locint' call, ensure that id being searched for is less
C than maximum (most times it won't be, saves searching entire list).
C
C Modified the node matching code to also index the nodes within the
C overlap region if doing a nodeset match. This has greatest benefit if
C not all nodes in the nodeset will be matched.
C
C Minor change in offset f -- 'dimension' to 'real'
C
C Revision 1.2  1992/09/03 20:35:02  gdsjaar
C Fixed up handling of mirror and offset - added a scale command which
C is combined with mirror, subroutine offset now does both offset and
C scale.  Added help to irennp.
C
c Revision 1.1  1992/09/02  22:58:05  gdsjaar
c Added mirroring and offsetting capability to gjoin.  Implemented at
c the equivalencing prompt.  Only affects the second database.
c
C=======================================================================
      subroutine offset (off, scale, crd, length)
      real crd(length)

      do 10 i=1, length
         crd(i) = scale * crd(i) + off
 10   continue
      
      return
      end
