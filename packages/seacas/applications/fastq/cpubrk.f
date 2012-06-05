C $Id: cpubrk.f,v 1.2 1992/09/09 16:50:49 gdsjaar Exp $
C $Log: cpubrk.f,v $
C Revision 1.2  1992/09/09 16:50:49  gdsjaar
C Added cpuifc function to cpubrk to remove need for suplib in fastq
C
c Revision 1.1.1.1  1990/11/30  11:05:27  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:05:25  gdsjaar
c Initial revision
c 
C
      LOGICAL FUNCTION CPUBRK(INPUT)
C***********************************************************************
C
C     FUNCTION CPUBRK = .TRUE. IF A CONTROL C HAS BEEN ENTERED AT TERMINAL
C
C***********************************************************************
      LOGICAL CPUIFC, INPUT
C
      IF (CPUIFC(INPUT)) THEN
         CPUBRK = .TRUE.
      ELSE
         CPUBRK = .FALSE.
      ENDIF
C
      RETURN
C
      END
C=======================================================================
      LOGICAL FUNCTION CPUIFC (LDUM)
C=======================================================================

C   --*** CPUIFC *** Dummy cancel function
C   --   Written by Amy Gilkey - revised 02/11/88
C   --
C   --CPUIFC returns the cancel flag as false.

      LOGICAL LDUM

      CPUIFC = .FALSE.

      RETURN
      END
