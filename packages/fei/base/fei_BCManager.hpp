#ifndef _fei_BCManager_hpp_
#define _fei_BCManager_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_defs.h"
#include "fei_macros.hpp"
#include "feiPoolAllocator.hpp"
#include <vector>

namespace fei {
  class Matrix;
}

class BCRecord;
class NodeDatabase;
class EqnBuffer;

/** Container for accumulating and organizing boundary-conditions.

    After BCManager is constructed, the 'addBCRecord' function is used to
    append boundary-condition records (each BCRecord is a specification of a
    boundary condition on a particular field on a particular node) to the
    BCManager's internal list of BCRecords. Once all boundary-conditions have
    been added, the 'consolidateBCs' function can be called to 'condense' the
    list of BCRecords, removing any duplicates.
    Note on removing duplicates:
    If duplicate dirichlet conditions are specified for a node/field pair, the
    last one added is kept and the others are discarded.
    If duplicate neumann conditions are specified for a node/field pair, the
    coefficients are added together.
 */

class BCManager {
 public:
   BCManager();
   virtual ~BCManager();

   void addBCRecords(int numNodes, const GlobalID* nodeIDs,
                    int fieldID, int fieldSize,
                    const double*const * alpha,
                    const double*const * beta,
                    const double*const * gamma);

   void addBCRecords(int idType, int numNodes, const GlobalID* nodeIDs,
                    int fieldID, int fieldSize,
                    const double*const * gamma,
                    const double*const * alpha);

   void addBCRecords(int idType, int numNodes, const GlobalID* nodeIDs,
                    int fieldID, int fieldSize,
                    const double*const * prescribedValues);

   void addBCRecords(int numNodes,
                    const GlobalID* nodeIDs,
                    int fieldID, int fieldSize,
                    const int* offsetsIntoField,
                    const double* prescribedValues);

   int finalizeBCEqns(fei::Matrix& matrix,
                      bool throw_if_bc_slave_conflict=false);

   int finalizeBCEqns(NodeDatabase& nodeDB,
                      EqnBuffer& bcEqns);

   int consolidateBCs();

   size_t getNumBCs();

   std::vector<const BCRecord*>& getBCRecords() { return( bcList_ ); }

   void clearAllBCs();

 private:
  BCManager(const BCManager& /*src*/);

  BCManager& operator=(const BCManager& /*src*/);

   std::vector<const BCRecord*> bcList_;
   feiPoolAllocator<BCRecord>* bcAlloc_;
   feiPoolAllocator<double>* coefAlloc_;
};

#endif

