C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C * Neither the name of NTESS nor the names of its
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

C=======================================================================
      SUBROUTINE NEWINI (IDNSUR, IDESUR, NSSUR, NUMATR)
C=======================================================================

C   $Id: newini.f,v 1.6 1999/02/17 15:26:56 gdsjaar Exp $
C   $Log: newini.f,v $
C   Revision 1.6  1999/02/17 15:26:56  gdsjaar
C   Converted to read/write exodusII database.
C   (copied repository from genshell to genshell2 prior to these changes)
C   Not tested much yet; probably still some problems
C
C   Revision 1.5  1991/04/10 22:26:53  gdsjaar
C   Fixed up node-set length problem
C
c Revision 1.4  1991/03/29  19:21:35  gdsjaar
c Removed the Gen3D variable blktyp
c
c Revision 1.3  1991/01/09  12:59:24  gdsjaar
c Initial conversion from GEN3D to GENSHELL, no BC yet
c
c Revision 1.2  90/10/01  15:39:50  gdsjaar
c Removed MAX() from dimension statement -- Non-ANSI usage
c 
c Revision 1.1.1.1  90/08/20  12:22:20  gdsjaar
c Gen3D Mesh Generation Program
c 
c Revision 1.1  90/08/20  12:22:19  gdsjaar
c Initial revision
c 

C   --*** NEWINI *** (GEN3D) Calculate 3D initial variables
C   --   Written by Amy Gilkey - revised 09/02/87
C   --
C   --NEWINI calculates the initial variables for the 3D database.
C   --The output number of nodes and elements and the length of the node
C   --sets and the side sets must be calculated before NEWINI is called.
C   --
C   --Parameters:
C   --   IDNSUR - IN - the number of surface node sets
C   --   IDESUR - IN - the number of surface side sets
C   --   NSSUR - IN - the number of nodes in the surface side set
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses NUMNP3, NUMEL3, LNPSN3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses LNPSNO, LESSEO, LESSNO of /DBNUM3/
C   --   Sets NUMNP3, NDIM3, NUMEL3, NELBL3,
C   --      NNPS3, LNPSN3, NESS3, LESSE3, LESSN3 of /DBNUM3/
C   --   Uses NNREPL, NEREPL of /PARAMS/

      include 'exodusII.inc'
      INCLUDE 'gs_dbtitl.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      INTEGER NUMATR(NELBLK)

C   --Database title - unchanged

      CONTINUE

C   --Number of dimensions

      NDIM3 = 3

C   --Number of nodes and elements - unchanged

      NUMNP3 = NUMNP
      NUMEL3 = NUMEL

C   --Number of element blocks

      NELBL3 = NELBLK

C   --Lengths of node sets set by NEWNPS
C   --Lengths of side sets set by NEWESS

C   --Number and lengths of sets, including front and back sets

      NNPS3 = NUMNPS + IDNSUR
      LNPSN3 = LNPSNO + IDNSUR*NUMNP
      NESS3 = NUMESS + IDESUR
      LESSE3 = LESSEO + IDESUR*NUMEL
      LESSN3 = LESSNO + IDESUR*NSSUR

C   --Number of attributes per block = 1

      DO 10 IBLK = 1, NELBLK
         NUMATR(IBLK) = 1
 10   CONTINUE

      RETURN
      END
