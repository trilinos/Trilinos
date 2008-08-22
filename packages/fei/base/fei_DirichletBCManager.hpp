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

#include <vector>

namespace fei {
class Matrix;

class DirichletBCManager {
 public:
  DirichletBCManager(){}
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

  size_t getNumBCRecords() const;

 private:
  std::vector<DirichletBCRecord> bcs_;
};//class DirichletBCManager
}//namespace fei
#endif

