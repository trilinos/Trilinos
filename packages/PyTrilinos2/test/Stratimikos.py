#!/usr/bin/env python
# @HEADER
# *****************************************************************************
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#
# Copyright 2022 NTESS and the PyTrilinos2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

from mpi4py import MPI
from argparse import ArgumentParser
try:
    # MPI, Timers, ParameterList
    from PyTrilinos2 import Teuchos
    # Linear algebra
    from PyTrilinos2 import Tpetra
    # Linear algebra interfaces used by Stratimikos
    from PyTrilinos2 import Thyra
    # Unified solver & preconditioner interface
    from PyTrilinos2 import Stratimikos
except ImportError:
    print("\nFailed to import PyTrilinos2. Consider setting the Python load path in your environment with\n export PYTHONPATH=${TRILINOS_BUILD_DIR}/packages/PyTrilinos2:${PYTHONPATH}\nwhere TRILINOS_BUILD_DIR is the build directory of your Trilinos build.\n")
    raise

try:
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams.update(mpl.rcParamsDefault)
    import matplotlib.pyplot as plt
    display = True
except:
    display = False


def assemble1DLaplacian(globalNumRows, comm):
    Tpetra_Map = Tpetra.Map()
    Tpetra_CrsGraph = Tpetra.CrsGraph()
    Tpetra_CrsMatrix = Tpetra.CrsMatrix()

    rowmap = Tpetra_Map(globalNumRows, 0, comm)

    graph = Tpetra_CrsGraph(rowmap, 3)
    for lid in range(rowmap.getMinLocalIndex(), rowmap.getMaxLocalIndex()+1):
        gid = rowmap.getGlobalElement(lid)
        indices = [gid]
        if gid > 0:
            indices.append(gid-1)
        if gid < rowmap.getMaxAllGlobalIndex():
            indices.append(gid+1)
        graph.insertGlobalIndices(gid, indices)
    graph.fillComplete()

    A = Tpetra_CrsMatrix(graph)
    for lid in range(rowmap.getMinLocalIndex(), rowmap.getMaxLocalIndex()+1):
        gid = rowmap.getGlobalElement(lid)
        indices = [gid]
        vals = [2.]
        if gid > 0:
            indices.append(gid-1)
            vals.append(-1.)
        if gid < rowmap.getMaxAllGlobalIndex():
            indices.append(gid+1)
            vals.append(-1.)
        A.replaceGlobalValues(gid, indices, vals)
    A.fillComplete()
    return A


def main():
    commMpi4py = MPI.COMM_WORLD
    comm = Teuchos.getTeuchosComm(commMpi4py)

    parser = ArgumentParser()
    parser.add_argument("--problemSize", default=1000, help="Global problem size", type=int)
    parser.add_argument("--solver", default="GMRES", choices=["LU", "CG", "BiCGStab", "GMRES"], help="Linear solver")
    parser.add_argument("--prec", default="None", choices=["None", "Jacobi", "Chebyshev", "ILU", "multigrid"], help="Preconditioner")
    parser.add_argument("--maxIts", default=100, help="Maximum number of iterations", type=int)

    # do a little dance to avoid redundant output
    doTerminate = False
    if comm.getRank() == 0:
        try:
            args = parser.parse_args()
        except SystemExit:
            doTerminate = True
        doTerminate = commMpi4py.bcast(doTerminate, root=0)
        if doTerminate:
            exit()
    else:
        doTerminate = commMpi4py.bcast(doTerminate, root=0)
        if doTerminate:
            exit()
        args = parser.parse_args()

    # Set up timers
    timer = Teuchos.StackedTimer("Main")
    Teuchos.TimeMonitor.setStackedTimer(timer)

    # Set up output streams
    cout = Teuchos.getCout()
    out = Teuchos.fancyOStream(cout)

    # Set up Tpetra matrices and vectors
    Tpetra_Map = Tpetra.Map()
    Tpetra_Vector = Tpetra.Vector()
    Tpetra_Export = Tpetra.Export()

    timer.start("Build matrix and vectors")
    tpetra_A = assemble1DLaplacian(args.problemSize, comm)
    tpetra_map = tpetra_A.getRowMap()

    n0 = (args.problemSize if comm.getRank() == 0 else 0)
    tpetra_map0 = Tpetra_Map(args.problemSize, n0, 0, comm)

    tpetra_x = Tpetra_Vector(tpetra_map, True)
    tpetra_b = Tpetra_Vector(tpetra_map, False)
    tpetra_residual = Tpetra_Vector(tpetra_map, False)

    tpetra_b.putScalar(1.)

    tpetra_A.describe(out)
    timer.stop("Build matrix and vectors")

    # Wrap Tpetra objects as Thyra objects
    thyra_map = Thyra.tpetraVectorSpace(tpetra_map)
    thyra_x = Thyra.tpetraVector(thyra_map, tpetra_x)
    thyra_b = Thyra.tpetraVector(thyra_map, tpetra_b)
    thyra_residual = Thyra.tpetraVector(thyra_map, tpetra_residual)
    thyra_A = Thyra.tpetraLinearOp(thyra_map, thyra_map, tpetra_A)

    # Set up linear solver
    linearSolverBuilder = Stratimikos.LinearSolverBuilder['double']()

    # Hook up preconditioners that are not enabled by default
    if hasattr(Stratimikos, "enableMueLu"):
        Stratimikos.enableMueLu(linearSolverBuilder, "MueLu")
    else:
        assert args.prec != "multigrid", "\"multigrid\" preconditioner requires the MueLu package."
    if hasattr(Stratimikos, "enableMueLuRefMaxwell"):
        Stratimikos.enableMueLuRefMaxwell(linearSolverBuilder, "MueLuRefMaxwell")
    if hasattr(Stratimikos, "enableMueLuMaxwell1"):
        Stratimikos.enableMueLuMaxwell1(linearSolverBuilder, "MueLuMaxwell1")

    # Print default parameters
    validParams = linearSolverBuilder.getValidParameters()
    if comm.getRank() == 0:
        print(validParams)

    params = Teuchos.ParameterList()

    # Set parameters for solver
    if args.solver == "LU":
        params["Linear Solver Type"] = "Amesos2"
    elif args.solver in ("CG", "GMRES", "BiCGStab"):
        if args.solver == "CG":
            BelosSolver = "Pseudo Block CG"
        elif args.solver == "GMRES":
            BelosSolver = "Block GMRES"
        elif args.solver == "BiCGStab":
            BelosSolver = "BiCGStab"
        params["Linear Solver Type"] = "Belos"
        params["Linear Solver Types"] = Teuchos.ParameterList()
        params["Linear Solver Types"]["Belos"] = Teuchos.ParameterList()
        params["Linear Solver Types"]["Belos"]["Solver Type"] = BelosSolver
        params["Linear Solver Types"]["Belos"]["Solver Types"] = Teuchos.ParameterList()
        params["Linear Solver Types"]["Belos"]["Solver Types"][BelosSolver] = Teuchos.ParameterList()
        params["Linear Solver Types"]["Belos"]["Solver Types"][BelosSolver]["Convergence Tolerance"] = 1e-8
        params["Linear Solver Types"]["Belos"]["Solver Types"][BelosSolver]["Maximum Iterations"] = args.maxIts
        params["Linear Solver Types"]["Belos"]["Solver Types"][BelosSolver]["Verbosity"] = 41
        params["Linear Solver Types"]["Belos"]["Solver Types"][BelosSolver]["Output Frequency"] = 1
        params["Linear Solver Types"]["Belos"]["Solver Types"][BelosSolver]["Output Style"] = 1

        params["Linear Solver Types"]["Belos"]["VerboseObject"] = Teuchos.ParameterList()
        params["Linear Solver Types"]["Belos"]["VerboseObject"]["Verbosity Level"] = "low"
    else:
        raise NotImplementedError("Unknown solver: {}".format(args.solver))

    # Set parameters for preconditioner
    if args.prec == "None":
        params["Preconditioner Type"] = "None"
    elif args.prec in ("Jacobi", "ILU", "Chebyshev"):
        params["Preconditioner Type"] = "Ifpack2"
        params["Preconditioner Types"] = Teuchos.ParameterList()
        params["Preconditioner Types"]["Ifpack2"] = Teuchos.ParameterList()
        if args.prec == "Jacobi":
            params["Preconditioner Types"]["Ifpack2"]["Prec Type"] = "relaxation"
            params["Preconditioner Types"]["Ifpack2"]["Ifpack2 Settings"] = Teuchos.ParameterList()
            params["Preconditioner Types"]["Ifpack2"]["Ifpack2 Settings"]["relaxation: type"] = "Jacobi"
            params["Preconditioner Types"]["Ifpack2"]["Ifpack2 Settings"]["relaxation: sweeps"] = 1
        elif args.prec == "Chebyshev":
            params["Preconditioner Types"]["Ifpack2"]["Prec Type"] = "chebyshev"
            params["Preconditioner Types"]["Ifpack2"]["Ifpack2 Settings"] = Teuchos.ParameterList()
            params["Preconditioner Types"]["Ifpack2"]["Ifpack2 Settings"]["chebyshev: degree"] = 2
        elif args.prec == "ILU":
            params["Preconditioner Types"]["Ifpack2"]["Prec Type"] = "ILUT"
    elif args.prec == "multigrid":
        params["Preconditioner Type"] = "MueLu"
    else:
        raise NotImplementedError("Unknown preconditioner: {}".format(args.prec))

    linearSolverBuilder.setParameterList(params)

    solverName = linearSolverBuilder.getLinearSolveStrategyName()
    precName = linearSolverBuilder.getPreconditionerStrategyName()
    solverFactory = Thyra.createLinearSolveStrategy(linearSolverBuilder)

    timer.start("Setup solver")
    thyra_invA = Thyra.linearOpWithSolve(solverFactory, thyra_A)
    timer.stop("Setup solver")
    assert thyra_invA.solveSupports(Thyra.NOTRANS)

    # Solve
    timer.start("Solve")
    status = thyra_invA.solve(Thyra.NOTRANS, thyra_b, thyra_x)
    timer.stop("Solve")

    # Compute residual
    tpetra_A.apply(tpetra_x, tpetra_residual)
    tpetra_residual.update(1, tpetra_b, -1)
    resNorm = tpetra_residual.norm2()

    if comm.getRank() == 0:
        print("Solver choice:                   ", args.solver)
        print("Stratimikos solver name:         ", solverName)
        print("Preconditioner choice:           ", args.prec)
        print("Stratimikos preconditioner name: ", precName)
    if solverName == "Belos":
        its = status.extraParameters["Iteration Count"]
        if comm.getRank() == 0:
            print('Norm of residual after {} iterations = {} '.format(its, resNorm))
    elif comm.getRank() == 0:
        print('Norm of residual = {} '.format(resNorm))

    # Gather solution vector on rank 0
    tpetra_x0 = Tpetra_Vector(tpetra_map0, True)
    export = Tpetra_Export(tpetra_map0, tpetra_map)
    tpetra_x0.doImport(source=tpetra_x, exporter=export, CM=Tpetra.CombineMode.REPLACE)

    # Print timings
    timer.stop("Main")
    options = Teuchos.StackedTimer.OutputOptions()
    options.output_fraction = options.output_histogram = options.output_minmax = True;
    timer.report(out, comm, options)

    if comm.getRank() == 0 and display:
        x0_view = tpetra_x0.getLocalViewHost()
        plt.figure()
        plt.plot(x0_view)
        plt.savefig('x0_view.png', dpi=800, bbox_inches='tight', pad_inches=0)

    comm.barrier()
    success = ((status.solveStatus == Thyra.ESolveStatus.SOLVE_STATUS_CONVERGED) and (resNorm < 1e-6))
    if comm.getRank() == 0:
        if success:
            print("OK")
        else:
            print("FAIL")


if __name__ == "__main__":
    # initialize kokkos
    defaultNode = Tpetra.Map.defaults['Node']
    if defaultNode in ('cuda', 'cuda_uvm', 'hip', 'hip_managed'):
        Tpetra.initialize_Kokkos(device_id=rank)
    else:
        Tpetra.initialize_Kokkos(num_threads=12)
    # Use a seperate function to make sure that all Tpetra objects get deallocated before we finalize Kokkos.
    main()
    # finalize kokkos
    Tpetra.finalize_Kokkos()
