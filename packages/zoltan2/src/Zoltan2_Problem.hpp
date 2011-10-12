#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Zoltan2_Model.hpp>

#include <Teuchos_CommHelpers.hpp>

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
  Problem(InputAdapter<User> &,
          Teuchos::ParameterList &params,
          const Teuchos::Comm<int> &comm = 
                       *(Teuchos::DefaultComm<int>::getComm().getRawPtr()));

  // Destructor
  virtual ~Problem() {};

  // Other methods
  virtual void solve() = 0;
  virtual void redistribute() = 0;

protected:
  RCP<InputAdapter<User> > inputAdapter_;
  RCP<Model<InputAdapter<User> > > model_;  
  RCP<Teuchos::ParameterList> params_;
  RCP<const Teuchos::Comm<int> > comm_;
  RCP<const Environment> env_;

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
  Teuchos::ParameterList &params,
  const Teuchos::Comm<int> &comm
) :
  inputAdapter_(rcpFromRef(input)),
  params_(rcpFromRef(params)),
  comm_(rcpFromRef(comm)),
  env_(Teuchos::RCP<Environment>(new Environment(params, comm_)))
{
  HELLO;
  cout << "KDDKDD input adapter type " << inputAdapter_->inputAdapterType() 
       << " " << inputAdapter_->inputAdapterName() 
       << " sizeof(scalar_t)= " 
       << sizeof(typename InputAdapter<User>::scalar_t) 
       << endl;
}

}

#endif
