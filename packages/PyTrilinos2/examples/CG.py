# @HEADER
# *****************************************************************************
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#
# Copyright 2022 NTESS and the PyTrilinos2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import unittest
from mpi4py import MPI

import numpy as np
from PyTrilinos2.PyTrilinos2 import Teuchos
from PyTrilinos2.PyTrilinos2 import Tpetra
from PyTrilinos2.PyTrilinos2 import MueLu
from PyTrilinos2.getTpetraTypeName import *
from math import sqrt

try:
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams.update(mpl.rcParamsDefault)
    import matplotlib.pyplot as plt
    display = True
except:
    display=False

def CG(A, x, b, max_iter=20, tol=1e-8, prec=None):
    r = type(b)(b, Teuchos.DataAccess.Copy)
    A.apply(x,r,Teuchos.ETransp.NO_TRANS,alpha=-1,beta=1)

    p = type(r)(r, Teuchos.DataAccess.Copy)
    q = type(r)(r, Teuchos.DataAccess.Copy)

    if prec is None:
        gamma = r.norm2()
    else:
        Br = type(r)(r, Teuchos.DataAccess.Copy)
        prec.apply(r, p)
        gamma = sqrt(r.dot(p))

    if gamma < tol:
        return 0
    for j in range(max_iter):
        A.apply(p, q)
        c = q.dot(p)
        alpha = gamma**2 / c
        x.update(alpha, p, 1)
        r.update(-alpha, q, 1)
        if prec is None:
            gamma_next = r.norm2()
            beta = gamma_next**2/gamma**2
            gamma = gamma_next
            if gamma < tol:
                return j+1
            p.update(1, r, beta)
        else:
            prec.apply(r, Br)
            gamma_next = sqrt(Br.dot(r))
            beta = gamma_next**2/gamma**2
            gamma = gamma_next
            if gamma < tol:
                return j+1
            p.update(1, Br, beta)
    return max_iter

def assemble1DLaplacian(n, comm):
    mapType = getTypeName('Map')
    graphType = getTypeName('CrsGraph')
    matrixType = getTypeName('CrsMatrix')

    mapT=mapType(n, 0, comm)
    graph = graphType(mapT, 3)
    for i in range(mapT.getMinLocalIndex(), mapT.getMaxLocalIndex()+1):
        global_i = mapT.getGlobalElement(i)
        indices = [global_i]
        if global_i > 0:
            indices.append(global_i-1)
        if global_i < mapT.getMaxAllGlobalIndex():
            indices.append(global_i+1)
        graph.insertGlobalIndices(global_i, indices)
    graph.fillComplete()

    A = matrixType(graph)
    for i in range(mapT.getMinLocalIndex(), mapT.getMaxLocalIndex()+1):
        global_i = mapT.getGlobalElement(i)
        indices = [global_i]
        vals = [2.]
        if global_i > 0:
            indices.append(global_i-1)
            vals.append(-1.)
        if global_i < mapT.getMaxAllGlobalIndex():
            indices.append(global_i+1)
            vals.append(-1.)
        A.replaceGlobalValues(global_i, indices, vals)
    A.fillComplete()

    return A


def main():
    comm = Teuchos.getTeuchosComm(MPI.COMM_WORLD)
    rank = comm.getRank()

    vectorType = getTypeName('Vector')

    n = 300000

    A = assemble1DLaplacian(n, comm)
    mapT = A.getRowMap()

    n0 = 0
    if rank == 0:
        n0 = n
    mapT0=type(mapT)(n, n0, 0, comm)

    x = vectorType(mapT, True)
    b = vectorType(mapT, False)
    residual = vectorType(mapT, False)

    b.putScalar(1.)

    p = Teuchos.ParameterList()
    P = MueLu.CreateTpetraPreconditioner(A, p)

    x.putScalar(0.)
    norm_x = x.norm2()
    if rank == 0:
        print('Norm of x before CG = {}'.format(norm_x))
    its = CG(A, x, b, max_iter=30, prec=P)
    norm_x = x.norm2()
    if rank == 0:
        print('Norm of x after {} iterations of CG = {} '.format(its, norm_x))

    A.apply(x, residual)
    residual.update(1, b, -1)
    resNorm = residual.norm2()
    if rank == 0:
        print('Norm of residual after {} iterations of CG = {} '.format(its, resNorm))

    x0 = vectorType(mapT0, True)
    export = getTypeName('Export')(mapT0, mapT)
    x0.doImport(source=x, exporter=export, CM=Tpetra.CombineMode.REPLACE)

    if rank == 0 and display:
        x0_view = x0.getLocalViewHost()
        plt.figure()
        plt.plot(x0_view)
        plt.savefig('x0_view.png', dpi=800, bbox_inches='tight',pad_inches = 0)

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if getDefaultNodeType() == 'cuda':
        Tpetra.initialize_Kokkos(device_id=rank)
    else:
        Tpetra.initialize_Kokkos(num_threads=12)
    main()
    Tpetra.finalize_Kokkos()
