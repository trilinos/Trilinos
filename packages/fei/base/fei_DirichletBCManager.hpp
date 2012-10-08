#ifndef _fei_DirichletBCManager_hpp_
#define _fei_DirichletBCManager_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_DirichletBCRecord.hpp>
#include <SNL_FEI_Structure.hpp>
#include <fei_VectorSpace.hpp>

#include <fei_Pool_alloc.hpp>
#include <map>

class NodeDatabase;
class EqnBuffer;

namespace fei {
class Matrix;

class DirichletBCManager {
 public:
  DirichletBCManager(SNL_FEI_Structure* structure)
   : structure_(structure), vecSpace_() {}

  DirichletBCManager(fei::SharedPtr<fei::VectorSpace> vecspace)
   : structure_(NULL), vecSpace_(vecspace) {}

  ~DirichletBCManager(){}

  void addBCRecords(int numBCs,
                    int IDType,
                    int fieldID,
                    int offsetIntoField,
                    const int* IDs,
                    const double* prescribedValues);

  void addBCRecords(int numBCs,
                    int IDType,
                    int fieldID,
                    const int* IDs,
                    const int* offsetsIntoField,
                    const double* prescribedValues);

  int finalizeBCEqns(fei::Matrix& matrix,
                     bool throw_if_bc_slave_conflict=false);

  int finalizeBCEqns(EqnBuffer& bcEqns);

  size_t getNumBCRecords() const;

  void clearAllBCs();

 private:
  int getEqnNumber(int IDType, int ID, int fieldID, int offsetIntoField);

  SNL_FEI_Structure* structure_;
  fei::SharedPtr<fei::VectorSpace> vecSpace_;

  typedef std::map<int,double,std::less<int>,
            fei_Pool_alloc<std::pair<const int, double> > > bc_map;
  bc_map bcs_;
};//class DirichletBCManager
}//namespace fei
#endif

