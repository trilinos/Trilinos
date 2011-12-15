#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>

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
template<typename Adapter>
class Problem {
public:
  
  // Constructors (there will be several to support novice interface)
  // Each will make sure the InputAdapter, parameters, etc. are set 
  // correctly before calling a common problem construction function.
  //KDDKDD How does simple interface work with Adapter template? Problem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &);
  //KDDKDD How does simple interface work with Adapter template? Problem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &, Teuchos::ParameterList &);
  Problem(Adapter *,
          Teuchos::ParameterList *params,
          const RCP<const Teuchos::Comm<int> > &comm = 
                       Teuchos::DefaultComm<int>::getComm());

  // Destructor
  virtual ~Problem() {};

  // Other methods
  virtual void solve() = 0;

protected:
  typedef typename Adapter::base_adapter_t base_adapter_t;

  RCP<const Adapter> inputAdapter_;
  RCP<const base_adapter_t> baseInputAdapter_;

  RCP<GraphModel<base_adapter_t> > graphModel_;  
  RCP<IdentifierModel<base_adapter_t> > identifierModel_;  

  RCP<const Model<base_adapter_t> > generalModel_;  

  // KDDKDD May want other models, too, for eval, printing, etc.
  RCP<Teuchos::ParameterList> params_;
  RCP<const Teuchos::Comm<int> > comm_;
  RCP<Environment> env_;
  RCP<const Environment> envConst_;

private:

};


#if 0 // KDDKDD How does simple interface work with Adapter template??
////////////////////////////////////////////////////////////////////////
//! Problem class constructor:  Tpetra matrix input must be converted
//! to XpetraMatrixAdapter.
template <typename Adapter>
Problem<Adapter>::Problem(
  Tpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> &A,
  Teuchos::ParameterList &p
) 
{
  HELLO;
  inputAdapter_ = rcp(new XpetraCrsMatrixInput<Z2PARAM_TEMPLATE>
                                (rcpFromRef(A)));
  params_ = rcpFromRef(p);
  cout << "KDDKDD input adapter type " << inputAdapter_->inputAdapterType() << " " << inputAdapter_->inputAdapterName() << endl;
}
#endif

template <typename Adapter>
Problem<Adapter>::Problem(
  Adapter *input,
  Teuchos::ParameterList *params,
  const RCP<const Teuchos::Comm<int> > &comm
) :
  inputAdapter_(RCP<Adapter>(input,false)),
  params_(RCP<Teuchos::ParameterList>(params,false)),
  comm_(comm),
  env_(Teuchos::RCP<Environment>(new Environment(*params, comm_)))
{
  HELLO;

  envConst_ = rcp_const_cast<const Environment>(env_);
  
//  cout << "KDDKDD input adapter type " << inputAdapter_->inputAdapterType() 
//       << " " << inputAdapter_->inputAdapterName() 
//       << " sizeof(scalar_t)= " 
//       << sizeof(typename Adapter::scalar_t) 
//       << endl;
}

}

#endif
