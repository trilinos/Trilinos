#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_config.h>
// TODO: SR We should have an input adapter, not sure why it was removed
//#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>
#include <Zoltan2_TemplateMacros.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Problem.hpp
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.

using Teuchos::RCP;
using Teuchos::rcp;
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
  // TODO: SR We should have an input adapter, not sure why it was removed
  //Problem(InputAdapter &);

  // Destructor
  virtual ~Problem() {};

  // Other methods
  virtual void solve() = 0;

protected:
  // TODO: SR We should have an input adapter, not sure why it was removed
  //RCP<InputAdapter> _inputAdapter;
  RCP<TpetraCrsMatrixInput<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > _inputAdapter;


  //TODO RCP<Model> _model;
  RCP<Teuchos::ParameterList> _params;

private:

};


////////////////////////////////////////////////////////////////////////
//! Problem class constructor:  Tpetra matrix input must be converted
//! to XpetraMatrixAdapter.
CONSISTENT_FN_TEMPLATE_LINE
Problem<CONSISTENT_TEMPLATE_PARAMS>::Problem(
  Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &A,
  Teuchos::ParameterList &p
) 
{
  HELLO;
  _inputAdapter = rcp(new TpetraCrsMatrixInput<CONSISTENT_TRILINOS_TEMPLATE_PARAMS>
                                (rcpFromRef(A)));
  _params = rcpFromRef(p);
}

}
#endif
