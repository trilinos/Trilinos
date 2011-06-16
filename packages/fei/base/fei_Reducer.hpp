/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Reducer_hpp_
#define _fei_Reducer_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_mpi.h>
#include <fei_Logger.hpp>
#include <fei_FillableVec.hpp>
#include <fei_FillableMat.hpp>
#include <fei_CSVec.hpp>
#include <fei_CSRMat.hpp>

namespace fei {
  class MatrixGraph;
  class Graph;
  class Matrix;
  class Vector;

  class Reducer : private fei::Logger {
   public:
    //@{ \name Constructors
    Reducer(fei::SharedPtr<FillableMat> globalSlaveDependencyMatrix,
            fei::SharedPtr<CSVec> g_vector,
            MPI_Comm comm);

    Reducer(fei::SharedPtr<fei::MatrixGraph> matrixGraph);
    //@}

    //@{ \name Destructor

    /** Destructor. */
    virtual ~Reducer();

    //@}
    void setLocalUnreducedEqns(const std::vector<int>& localUnreducedEqns);

    
    /** Set the matrix-graph structure. This is the nonzero structure for
        locally-owned matrix rows.
    */
    void addGraphEntries(fei::SharedPtr<fei::SparseRowGraph> matrixGraph);

    void addGraphIndices(int numRows, const int* rows,
                         int numCols, const int* cols,
                         fei::Graph& graph);

    void addSymmetricGraphIndices(int numIndices, const int* indices,
                                  bool diagonal,
                                  fei::Graph& graph);

    /** Put a C-style table (array of pointers) of coefficient data into the
        matrix.  This is a rectangular array of coefficients for
        rows/columns defined by the 'rows' and 'cols' lists.
        If the sum_into argument is true, values should be added to any that
        already exist at the specified locations. Otherwise (if sum_into is
        false) incoming values should overwrite already-existing values.
     */
    int addMatrixValues(int numRows, const int* rows,
                        int numCols, const int* cols,
                        const double* const* values,
                        bool sum_into,
                        fei::Matrix& feimat,
                        int format);

    /** Put coefficient data into a vector at the specified global indices.
      If any specified indices are out of range (negative or too large) the
      corresponding positions in the values array will not be referenced,
      and a positive warning code will be returned.

      @param numValues Length of caller-allocated 'globalIndices' and
                       'values' arrays.

      @param globalIndices List of global-indices specifying the locations in
                the vector for incoming values to be placed.

      @param values List of incoming values.

      @param sum_into If true, incoming values should be added to values that
                 may already be in the specified locations. If sum_into is
                 false, then incoming values should overwrite existing values.

      @param soln_vector If true, incoming values should be placed in the
                   solution vector. Otherwise, they should be placed in the
                   rhs vector.

      @param vectorIndex If the linear system has multiple rhs/soln vectors,
                       then this parameter specifies which vector the incoming
                       values should be put into.
    */
    int addVectorValues(int numValues,
                        const int* globalIndices,
                        const double* values,
                        bool sum_into,
                        bool soln_vector,
                        int vectorIndex,
                        fei::Vector& feivec);

    int copyOutVectorValues(int numValues,
                             const int* globalIndices,
                             double* values,
                             bool soln_vector,
                             int vectorIndex,
                             fei::Vector& feivec);

    void getSlaveMasterEqns(int slaveEqn, std::vector<int>& masterEqns);
    bool isSlaveEqn(int unreducedEqn) const;
    bool isSlaveCol(int unreducedEqn) const;

    /** Given an equation-number in the caller's unreduced index-space,
      return the corresponding equation in the reduced space.
      If unreducedEqn is a slave, an exception will be thrown.
    */
    int translateToReducedEqn(int unreducedEqn) const;
    int translateFromReducedEqn(int reduced_eqn) const;
    void assembleReducedGraph(fei::Graph* graph,
                              bool global_gather=true);
    void assembleReducedGraph(fei::SparseRowGraph* srgraph);
    void assembleReducedMatrix(fei::Matrix& matrix);
    void assembleReducedVector(bool soln_vector,
                               fei::Vector& feivec);

    std::vector<int>& getLocalReducedEqns();

    void initialize();
   private:
    void expand_work_arrays(int size);

    fei::CSRMat csrD_;
    int* slavesPtr_;
    fei::FillableMat Kii_, Kid_, Kdi_, Kdd_;
    fei::CSRMat csrKii, csrKid, csrKdi, csrKdd;
    fei::FillableVec fi_, fd_;
    fei::CSVec csfi, csvec, csvec_i;
    fei::CSRMat tmpMat1_, tmpMat2_;
    fei::CSVec tmpVec1_, tmpVec2_;

    fei::CSVec csg_;
    bool g_nonzero_;

    std::vector<int> localUnreducedEqns_;
    std::vector<int> localReducedEqns_;
    std::vector<int> nonslaves_;
    std::vector<int> reverse_;
    bool* isSlaveEqn_;
    int numGlobalSlaves_;
    int numLocalSlaves_;
    int firstLocalReducedEqn_;
    int lastLocalReducedEqn_;
    int lowestGlobalSlaveEqn_;
    int highestGlobalSlaveEqn_;

    int localProc_;
    int numProcs_;
    MPI_Comm comm_;
    std::string dbgprefix_;
    unsigned mat_counter_;
    unsigned rhs_vec_counter_;

    bool* bool_array_;
    int* int_array_;
    double* double_array_;
    int array_len_;

    std::vector<double> work_1D_;
    std::vector<const double*> work_2D_;
  };//class Reducer

}//namespace fei

#endif // _fei_Reducer_hpp_

