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
from PyTrilinos2.getTpetraTypeName import getTypeName
from math import sqrt


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


class TestCG(unittest.TestCase):
    def test_all(self):
        comm = Teuchos.getTeuchosComm(MPI.COMM_WORLD)

        mapType = getTypeName('Map')
        graphType = getTypeName('CrsGraph')
        matrixType = getTypeName('CrsMatrix')
        vectorType = getTypeName('Vector')
        multivectorType = getTypeName('MultiVector')

        n = 3000

        mapT=mapType(n, 0, comm)
        print(mapT)
        print('mapT.getMinLocalIndex() = '+str(mapT.getMinLocalIndex()))
        print('mapT.getMaxLocalIndex() = '+str(mapT.getMaxLocalIndex()))
        print('mapT.getMinGlobalIndex() = '+str(mapT.getMinGlobalIndex()))
        print('mapT.getMaxGlobalIndex() = '+str(mapT.getMaxGlobalIndex()))
        print('Tpetra.getDefaultComm().getSize() = '+str(Tpetra.getDefaultComm().getSize()))
        mv=multivectorType(mapT, 3, True)
        mv.replaceLocalValue(0,0,1.23)
        mv.replaceLocalValue(0,1,1.23)
        #mv.randomize(0,-2)
        v0=mv.getVector(0)
        v1=mv.getVector(1)
        print(mv.description())
        print(v0.description())
        print(v0.norm2())
        print(v0.dot(v1))

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

        print(A.getGlobalNumEntries())
        print(A.description())
        print(A.getFrobeniusNorm())

        print(v0.norm2())
        print(v1.norm2())
        A.apply(v0,v1)
        print(v1.norm2())


        x = vectorType(mapT, True)
        b = vectorType(mapT, False)
        residual = vectorType(mapT, False)

        b.randomize(0,-2)

        print('Norm of x before CG = {}'.format(x.norm2()))
        print('Norm of b = '+str(b.norm2()))
        its = CG(A, x, b, max_iter=n)
        print('Norm of x after {} iterations of CG = {} '.format(its, x.norm2()))

        A.apply(x, residual)
        residual.update(1, b, -1)
        resNorm = residual.norm2()
        print('Norm of residual after {} iterations of CG = {} '.format(its, resNorm))

        self.assertAlmostEqual(resNorm, 0., delta=1e-5)

        p = Teuchos.ParameterList()
        P = MueLu.CreateTpetraPreconditioner(A, p)

        x.putScalar(0.)
        print('Norm of x before CG = {}'.format(x.norm2()))
        its = CG(A, x, b, max_iter=30, prec=P)
        print('Norm of x after {} iterations of CG = {} '.format(its, x.norm2()))

        A.apply(x, residual)
        residual.update(1, b, -1)
        resNorm = residual.norm2()
        print('Norm of residual after {} iterations of CG = {} '.format(its, resNorm))

        self.assertAlmostEqual(resNorm, 0., delta=1e-5)


if __name__ == '__main__':
    unittest.main()
