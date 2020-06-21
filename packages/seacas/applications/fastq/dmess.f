C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: dmess.f,v 1.2 2007/07/24 13:10:18 gdsjaar Exp $
C $Log: dmess.f,v $
C Revision 1.2  2007/07/24 13:10:18  gdsjaar
C Fix problem with boundary condition memory overwrite.
C
C Remove old ls5 and r25 terminal tests
C
C Revision 1.1.1.1  1990/11/30 11:06:17  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:06:16  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]DMESS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DMESS (DEV1, TEXT)
C***********************************************************************
C
C  SUBROUTINE DMESS = PRINTS A ONE LINE MESSAGE AT THE BOTTOM OF THE
C                       SCREEN
C
C***********************************************************************
C
      CHARACTER*(*) TEXT, DEV1*3
C
      CALL MESAGE (TEXT)
      RETURN
C
      END
