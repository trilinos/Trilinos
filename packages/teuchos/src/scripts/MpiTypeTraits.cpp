void safeMpiDatatypeFree (MPI_Datatype* datatype) {
  int finalized = 0;
  const int err = MPI_Finalized (&finalized);
  if (err == MPI_SUCCESS && ! finalized && *datatype != MPI_DATATYPE_NULL) {
    // Don't throw an exception on error, since this function will be
    // called in a destructor.  Just ignore any errors and move on.
    // If successful, MPI_Type_free sets its input to
    // MPI_DATATYPE_NULL.
    (void) MPI_Type_free (datatype);
  }
}
