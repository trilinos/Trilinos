// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __test_utils_hpp__
#define __test_utils_hpp__

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Teko_Utilities.hpp"

#ifdef TEKO_HAVE_EPETRA
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#endif
#include "Teuchos_Comm.hpp"

#include "Tpetra_Vector.hpp"

#include <iostream>
#include <list>

namespace Teko {
namespace Test {

// build a 2x2 matrix...only in serial
#ifdef TEKO_HAVE_EPETRA
const Teuchos::RCP<const Thyra::LinearOpBase<double> > build2x2(const Epetra_Comm& comm, double a,
                                                                double b, double c, double d);
#endif
const Teuchos::RCP<const Thyra::LinearOpBase<ST> > build2x2(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm, ST a, ST b, ST c, ST d);

// prints a vector, with string "s" as the name
void Print(std::ostream& os, const std::string& s,
           const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& v);

// prints a vector, with string "s" as the name
inline std::string Print(const std::string& s,
                         const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& v) {
  std::stringstream ss;
  Print(ss, s, v);
  return ss.str();
}

// builds a single vector of the identity matrix
Teuchos::RCP<Thyra::VectorBase<double> > BuildIVector(
    int j, const Teuchos::RCP<const Thyra::VectorSpaceBase<double> >& vs);

// peforms an extraction of the matrix, in addition to printing it
// DANGER! This is a CPU intense routine, should only be done on small, in-core matricies
void HardExtract(std::ostream& os, const Teuchos::RCP<const Thyra::LinearOpBase<double> >& A);

// add two Thyra vectors
const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(
    const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& x,
    const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& y);
const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(
    double ax, const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& x, double ay,
    const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& y);

// compute ||x-y||_2
double Difference(const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& x,
                  const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& y);

// construct a diagonal matrix
#ifdef TEKO_HAVE_EPETRA
const Teuchos::RCP<const Thyra::LinearOpBase<double> > DiagMatrix(int cnt, double* vec,
                                                                  std::string label = "");
#endif

const Teuchos::RCP<const Thyra::LinearOpBase<ST> > DiagMatrix_tpetra(GO cnt, ST* vec,
                                                                     std::string label = "");

// 2-Vector
#ifdef TEKO_HAVE_EPETRA
const Teuchos::RCP<const Thyra::MultiVectorBase<double> > BlockVector(
    const Epetra_Vector& u, const Epetra_Vector& v,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double> >& vs);
#endif

const Teuchos::RCP<const Thyra::MultiVectorBase<ST> > BlockVector(
    const Tpetra::Vector<ST, LO, GO, NT>& u, const Tpetra::Vector<ST, LO, GO, NT>& v,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<ST> >& vs);

class UnitTest {
 public:
  virtual ~UnitTest() {}
  virtual void initializeTest()      = 0;
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                      int& totalrun) = 0;
  virtual bool isParallel() const    = 0;

  static void AddTest(const Teuchos::RCP<UnitTest>& ut, const std::string& name);
#ifdef TEKO_HAVE_EPETRA
  static bool RunTests(int verbosity, std::ostream& stdstrm, std::ostream& failstrm);
#endif
  static bool RunTests_tpetra(int verbosity, std::ostream& stdstrm, std::ostream& failstrm);

#ifdef TEKO_HAVE_EPETRA
  static Teuchos::RCP<const Epetra_Comm> GetComm();
#endif
  static Teuchos::RCP<const Teuchos::Comm<int> > GetComm_tpetra();

#ifdef TEKO_HAVE_EPETRA
  static void SetComm(const Teuchos::RCP<const Epetra_Comm>& c);
#endif
  static void SetComm_tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >& c);

  static void ClearTests();

 protected:
  static std::list<std::pair<Teuchos::RCP<UnitTest>, std::string> > testList;
#ifdef TEKO_HAVE_EPETRA
  static Teuchos::RCP<const Epetra_Comm> comm_;
#endif
  static Teuchos::RCP<const Teuchos::Comm<int> > comm_tpetra_;

#ifdef TEKO_HAVE_EPETRA
  static bool CheckParallelBools(bool myBool, int& failPID);
#endif
  static bool CheckParallelBools_tpetra(bool myBool, int& failPID);
};

inline const std::string toString(bool status) { return status ? "PASSED" : "FAILED"; }

}  // namespace Test
}  // end namespace Teko

#define Teko_ADD_UNIT_TEST(str, name) Teko::Test::UnitTest::AddTest(Teuchos::rcp(new str()), #name)
#define Teko_TEST_MSG(os, level, msgp, msgf)                     \
  {                                                              \
    int failPID = -1;                                            \
    status      = UnitTest::CheckParallelBools(status, failPID); \
    if (verbosity >= level && status)                            \
      os << msgp << std::endl;                                   \
    else if (verbosity >= level && not status)                   \
      os << msgf << " ( PID = " << failPID << " )" << std::endl; \
  }

#define Teko_TEST_MSG_tpetra(os, level, msgp, msgf)                     \
  {                                                                     \
    int failPID = -1;                                                   \
    status      = UnitTest::CheckParallelBools_tpetra(status, failPID); \
    if (verbosity >= level && status)                                   \
      os << msgp << std::endl;                                          \
    else if (verbosity >= level && not status)                          \
      os << msgf << " ( PID = " << failPID << " )" << std::endl;        \
  }

#define TEST_EQUALITY(x, y, msg)       \
  status = (x == y);                   \
  if (not status || verbosity >= 10) { \
    os << msg << std::endl;            \
  }                                    \
  allPassed &= status;

#define TEST_NOT_EQUAL(x, y, msg)      \
  status = (x != y);                   \
  if (not status || verbosity >= 10) { \
    os << msg << std::endl;            \
  }                                    \
  allPassed &= status;

#define TEST_MSG(msg)       \
  if (verbosity >= 10) {    \
    os << msg << std::endl; \
  }

#define TEST_ASSERT(x, msg)            \
  status = x;                          \
  if (not status || verbosity >= 10) { \
    os << msg << std::endl;            \
  }                                    \
  allPassed &= status;

#endif
