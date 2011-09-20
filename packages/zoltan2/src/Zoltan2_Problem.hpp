#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_config.h>
#include <Zoltan2_InputAdapter.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>

#include <Kokkos_DefaultNode.hpp>

////////////////////////////////////////////////////////////////////////
#define HELLO cout << "Hello from " << __func__ << endl

#define CONSISTENT_CLASS_TEMPLATE_LINE \
        template <typename Scalar=float, \
                  typename LNO=int, typename GNO=int, typename LID=LNO, typename GID=GNO, \
                  typename Node=Kokkos::DefaultNode::DefaultNodeType>

#define CONSISTENT_FN_TEMPLATE_LINE \
        template <typename Scalar, \
                  typename LNO, typename GNO, typename LID, typename GID, \
                  typename Node>

#define CONSISTENT_TEMPLATE_PARAMS \
        Scalar, LNO, GNO, LID, GID, Node

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
  Problem(InputAdapter &);

  ~Problem() {};

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
  Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> *A, 
  Teuchos::ParameterList *p
) 
{
  HELLO;
  //TODO _inputAdapter = rcp(new TpetraMatrixAdapter(rcp(A))); 
  _params = rcp(p);
}

}
#endif
