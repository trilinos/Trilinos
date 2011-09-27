#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>
#include <Zoltan2_InputAdapter.hpp>

#include <Tpetra_CrsMatrix.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Problem.hpp
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.

using std::cout;
using std::endl;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \class Problem
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.
CONSISTENT_CLASS_TEMPLATE_LINE
class Problem {

public:
  // Constructors (there will be several to support novice interface)
  // Each will make sure the InputAdapter, parameters, etc. are set 
  // correctly before calling a common problem construction function.
  Problem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &);
  Problem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &, Teuchos::ParameterList &);
  Problem(InputAdapter &);

  // Destructor
  virtual ~Problem() {};

  // Other methods
  virtual void solve() = 0;

protected:
  RCP<InputAdapter> _inputAdapter;
  //TODO RCP<Model> _model;
  RCP<Teuchos::ParameterList> _params;

private:

};


////////////////////////////////////////////////////////////////////////
//! Problem class constructor:  Tpetra matrix input must be converted
//! to XpetraMatrixAdapter.
CONSISTENT_FN_TEMPLATE_LINE
Problem<CONSISTENT_TEMPLATE_PARAMS>::Problem(
  Tpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> &A,
  Teuchos::ParameterList &p
) 
{
  HELLO;
  _inputAdapter = rcp(new TpetraCrsMatrixInput<CONSISTENT_TEMPLATE_PARAMS>
                                (rcpFromRef(A)));
  _params = rcpFromRef(p);
}

}
#endif
