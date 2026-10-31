def printP0(*args, comm=None, root=0, **kwargs):
    """
    Print only from the root rank of an MPI communicator.

    If comm is not provided, attempt to use MPI.COMM_WORLD.
    If mpi4py is unavailable, behave like normal print().
    """
    if comm is None:
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
        except ImportError:
            print(*args, **kwargs)
            return

    if comm.Get_rank() == root:
        print(*args, **kwargs)
        
