//  MPI_Reduce_scatter( msg_count, &nrecvs_, counts, MPI_INT, MPI_SUM, comm_ );
  MPI_Reduce(msg_count, counts, nprocs, MPI_INT, MPI_SUM, 0, comm_);
  MPI_Scatter(counts, 1, MPI_INT, &nrecvs_, 1, MPI_INT, 0, comm_);
