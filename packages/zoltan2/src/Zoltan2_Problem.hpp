#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Zoltan2_Model.hpp>

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
template<Z2CLASS_TEMPLATE>
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
  virtual void redistribute() = 0;

protected:
  RCP<InputAdapter> inputAdapter_;
  //TODO RCP<Model<InputAdapter, Z2FN_TEMPLATE> > model_;  
  RCP<Teuchos::ParameterList> params_;

private:

};


////////////////////////////////////////////////////////////////////////
//! Problem class constructor:  Tpetra matrix input must be converted
//! to XpetraMatrixAdapter.
template <Z2FN_TEMPLATE>
Problem<Z2PARAM_TEMPLATE>::Problem(
  Tpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> &A,
  Teuchos::ParameterList &p
) 
{
  HELLO;
  inputAdapter_ = rcp(new TpetraCrsMatrixInput<Z2PARAM_TEMPLATE>
                                (rcpFromRef(A)));
  params_ = rcpFromRef(p);
  cout << "KDDKDD input adapter type " << inputAdapter_->inputAdapterType() << " " << inputAdapter_->inputAdapterName() << endl;
}

}
#endif
