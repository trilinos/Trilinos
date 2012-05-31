C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
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

      block data cdrabc
      integer kwrtfl,krdfl,koutfl,kinfl,kwrdsz,kbytel,kcpw,kbaud,
     1kcomtp
      common /cdrcom/ kwrtfl,krdfl,koutfl,kinfl,kwrdsz,kbytel,kcpw,
     1 kbaud,kcomtp
      common /cdrcm2 / KGNAME
      character*80 KGNAME
      integer kuntfd(0:999)
      common /cdrunx/ kuntfd
      logical onode
      common /vconod/ onode
      integer machin(3),maclen
      integer kidsiz,kjobid(4),kusrsz,kusrid(4),kszrou
      integer kjrout(4),ksecur,kjtime(3),kjdate(3)
      common /vcjob/ kidsiz,kjobid,kusrsz,kusrid,kszrou,
     1               kjrout,ksecur,kjtime,kjdate,machin,maclen

      data kwrtfl /6/
      data krdfl /5/
      data koutfl /77/
      data kinfl /66/
      data kwrdsz /64/
      data kbytel /8/
      data kcpw /8/
      data kbaud /9600/

c  COMPUTER TYPE WHERE DATA ORIGINATED  
c  VAX NETWORK = 4                     
c  NOS         = 5                    
c  SCOPE 2.1   = 6                   
c  UNIVAC      = 7                  
c  cray        = 8                 
c  ALLIANT     = 9                
c  APOLLO      = 10              
c  SUN         = 11             
c  HP          = 12            
c  DEC/ULTRIX  = 13

      data kcomtp /11/

c initialize unit to file descriptor table (cdrunx)
      data kuntfd/1000 * 0/

C initialize graphics output filename to blank (cdrcm2)
      data kgname/' '/

      end

