
      program simple

      include 'mpif.h'

      integer info
      integer np
      integer rank

      call mpi_init(info)
      call mpi_comm_size(mpi_comm_world,np,info)
      call mpi_comm_rank(mpi_comm_world,rank,info)
      print *,'rank ',rank,' of ',np,' initialized'
      call mpi_finalize(info)

      end
