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
template<typename User>
class Problem {
public:
  
  // Constructors (there will be several to support novice interface)
  // Each will make sure the InputAdapter, parameters, etc. are set 
  // correctly before calling a common problem construction function.
  //KDDKDD How does simple interface work with User template? Problem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &);
  //KDDKDD How does simple interface work with User template? Problem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &, Teuchos::ParameterList &);
  Problem(InputAdapter<User> &);

  // Destructor
  virtual ~Problem() {};

  // Other methods
  virtual void solve() = 0;
  virtual void redistribute() = 0;

protected:
  RCP<InputAdapter<User> > inputAdapter_;
  RCP<Model<InputAdapter<User> > > model_;  
  RCP<Teuchos::ParameterList> params_;

private:

};


#if 0 // KDDKDD How does simple interface work with User template??
////////////////////////////////////////////////////////////////////////
//! Problem class constructor:  Tpetra matrix input must be converted
//! to XpetraMatrixAdapter.
template <typename User>
Problem<User>::Problem(
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
#endif

template <typename User>
Problem<User>::Problem(
  InputAdapter<User> &input,
  Teuchos::ParameterList &params
)
{
  HELLO;
  inputAdapter_ = rcp(input);
  params_ = rcpFromRef(p);
  cout << "KDDKDD input adapter type " << inputAdapter_->inputAdapterType() 
       << " " << inputAdapter_->inputAdapterName() << endl;
}

}

}
#endif
