#ifndef _InputFileReader_hpp_
#define _InputFileReader_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_ParameterSet.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Factory.hpp>
#include "fei_Record.hpp"
#include <snl_fei_Constraint.hpp>

#include <test_utils/feitester.hpp>

class InputFileReader : public feitester {
 public:
  InputFileReader(MPI_Comm comm,
                  fei::SharedPtr<fei::Factory> factory);
  virtual ~InputFileReader();

  int readInputFile(const char* fileName);

  int readParameters(const char* fileName);

  void setSolnFileName(const std::string& solnFileName)
  { solnFileName_ = solnFileName; }

  void setCheckFileName(const std::string& checkFileName)
  { checkFileName_ = checkFileName; }

  //----feitester methods: --------
  const char* getName() { return( "InputFileReader"); }

  int testInitialization();

  int testLoading();

  int testSolve();

  int testCheckResult();

  void dumpMatrixFiles();

  void setParameter(const char* param);

  //----end of feitester methods-----

  fei::SharedPtr<fei::ParameterSet> feiParameterSet;

  fei::SharedPtr<fei::VectorSpace>& getRowSpace();

  fei::SharedPtr<fei::VectorSpace> feiVectorSpace_col;

  fei::SharedPtr<fei::MatrixGraph>& getMatrixGraph();

  fei::SharedPtr<fei::Matrix>& getMatrix();

  fei::SharedPtr<fei::Vector>& getVector(bool soln);

  fei::SharedPtr<fei::LinearSystem>& getLinearSystem();

 private:
  void readData(FEI_ISTREAM* instr, const char* keyword);

  void readField(FEI_ISTREAM* instr);
  void readIDType(FEI_ISTREAM* instr);
  void readSimplePattern(FEI_ISTREAM* instr);
  void readSymmElemBlock(FEI_ISTREAM* instr);
  void readNonsymmConnBlock(FEI_ISTREAM* instr);
  void readSymmConnectivity(FEI_ISTREAM* instr);
  void readNonsymmConnectivity(FEI_ISTREAM* instr);
  void readSlaveConstraint(FEI_ISTREAM* instr);
  void readElemStiffness(FEI_ISTREAM* instr);
  void readElemLoad(FEI_ISTREAM* instr);
  void readLagrangeConstraint(FEI_ISTREAM* instr);
  void readNonsymmCoefficients(FEI_ISTREAM* instr);
  void readSharedIDs(FEI_ISTREAM* instr);
  void loadLagrangeConstraints(fei::SharedPtr<fei::LinearSystem> linsys);
  void save_soln(int idType, const char* fileName);

  MPI_Comm comm_;

  fei::SharedPtr<fei::Factory> factory_;
  fei::SharedPtr<fei::VectorSpace> feiVectorSpace_row_;
  fei::SharedPtr<fei::MatrixGraph> feiMatrixGraph_;
  bool initCompleteCalled_;
  fei::SharedPtr<fei::Matrix> feiMatrix_;
  fei::SharedPtr<fei::Vector> feiVector_soln_;
  fei::SharedPtr<fei::Vector> feiVector_rhs_;
  fei::SharedPtr<fei::LinearSystem> feiLinearSystem_;

  std::string solnFileName_;
  std::string checkFileName_;
  int solveCounter_;

  typedef snl_fei::Constraint<fei::Record*> ConstraintType;

  std::map<int,ConstraintType*> lagrangeConstraints_;
  int crID_;
};

#endif

