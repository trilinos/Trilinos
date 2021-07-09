C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
C==================================================
C Parameters for the NE_TEST fortran interface test
C program
C==================================================

C Constants for init_global functions
      PARAMETER (NNG = 100)
      PARAMETER (NEG = 50)
      PARAMETER (NEBG = 5)
      PARAMETER (NNSG = 2)
      PARAMETER (NSSG = 3)

C Constants for load-balance functions
      PARAMETER (NPROC = 10)
      PARAMETER (NPROCF = NPROC)
      PARAMETER (NINTN = 200)
      PARAMETER (NBORN = 10)
      PARAMETER (NEXTN = 5)
      PARAMETER (NINTE = 100)
      PARAMETER (NBORE = 20)
      INTEGER NNCMAP, NECMAP
      PARAMETER (NNCMAP = 4)
      PARAMETER (NECMAP = 2)

C Constants for communication map
      PARAMETER (NCNTCM = 20)
      INTEGER ECNTCM
      PARAMETER (ECNTCM = 17)

C Some multiplied values
C   NNCMAP * NPROCF
      PARAMETER (NNCXNPF = 40)
C   NECMAP * NPROCF
      PARAMETER (NECXNPF = 20)
