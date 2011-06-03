/**
 * \file   Solver_Test.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Wed May 25 12:17:25 2011
 * 
 * \brief  Tests Amesos2 solver interfaces using various matrix/vector
 *         objects, scalar/ordinal types, and input matrices.
 */

/*
 * Use an xml file to specify:
 *
 * input matrix
 * - Solver to test with that matrix
 * - Matrix/Multivector objects to use (Epetra, Tpetra, ...)
 * - Solver parameters (Transpose?, etc..)
 */
