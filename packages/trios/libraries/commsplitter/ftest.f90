
      PROGRAM commsplitter_test

      INCLUDE 'mpif.h'

      INTEGER info
      INTEGER np
      INTEGER rank
      CHARACTER*(1024) hostname
      INTEGER hostname_len

      CALL MPI_Init(info)
      CALL MPI_Comm_size(mpi_comm_world,np,info)
      CALL MPI_Comm_rank(mpi_comm_world,rank,info)
      CALL MPI_Get_processor_name(hostname, hostname_len, info);
      PRINT *,'ftest: rank ',rank,' of ',np,' executing on ', trim(hostname)
      CALL MPI_Finalize(info)

      end
