// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_OrderingProblem.hpp
    \brief Defines the OrderingProblem class.
*/

#ifndef _ZOLTAN2_ORDERINGPROBLEM_HPP_
#define _ZOLTAN2_ORDERINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
//#include <Zoltan2_OrderingAlgorithms.hpp> // TODO: Fix include path?
#include "algorithms/order/Zoltan2_OrderingAlgorithms.hpp"
#include <Zoltan2_OrderingSolution.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <string>
#ifdef HAVE_ZOLTAN2_OVIS
#include <ovis.h>
#endif

using Teuchos::rcp_dynamic_cast;
using namespace std;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////

/*! \brief OrderingProblem sets up ordering problems for the user.
 *
 *  The OrderingProblem is the core of the Zoltan2 ordering API.
 *  Based on the the user's input and parameters, the OrderingProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the OrderingProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo follow ordering with partitioning
 */

template<typename Adapter>
class OrderingProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lno_t lno_t;

#ifdef HAVE_ZOLTAN2_MPI
   typedef Teuchos::OpaqueWrapper<MPI_Comm> > mpiWrapper_t;
#endif

  /*! \brief Destructor
   */
  virtual ~OrderingProblem() {};


#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor that takes an MPI communicator
   */
  OrderingProblem(Adapter *A, ParameterList *p, MPI_Comm comm) 
                      : Problem<Adapter>(A, p, comm) 
  {
    HELLO;
    createOrderingProblem();
  };
#endif

  /*! \brief Constructor that uses a default communicator
   */
  OrderingProblem(Adapter *A, ParameterList *p) : Problem<Adapter>(A, p) 
  {
    HELLO;
    createOrderingProblem();
  };

  //!  \brief Direct the problem to create a solution.
  //
  //    \param updateInputData   If true this indicates that either
  //          this is the first attempt at solution, or that we
  //          are computing a new solution and the input data has
  //          changed since the previous solution was computed.
  //          If false, this indicates that we are computing a
  //          new solution using the same input data was used for
  //          the previous solution, even though the parameters
  //          may have been changed.
  //
  //  For the sake of performance, we ask the caller to set \c updateInputData
  //  to false if he/she is computing a new solution using the same input data,
  //  but different problem parameters, than that which was used to compute
  //  the most recent solution.
  
  void solve(bool updateInputData=true);

  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  OrderingSolution<gid_t, lno_t> *getSolution() {
    return solution_.getRawPtr();
  };

private:
  void createOrderingProblem();

  RCP<OrderingSolution<gid_t, lno_t> > solution_;

  RCP<Comm<int> > problemComm_;

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpiComm_;
#endif
};

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void OrderingProblem<Adapter>::solve(bool newData)
{
  HELLO;

  size_t nVtx = this->graphModel_->getLocalNumVertices();

  // TODO: Assuming one MPI process now. nVtx = ngids = nlids
  this->solution_ = rcp(new OrderingSolution<gid_t, lno_t>(nVtx, nVtx));

  // Determine which algorithm to use based on defaults and parameters.
  // TODO: Use RCM if graph model is defined, otherwise use Natural.
  // Need some exception handling here, too.

  string method = this->params_->template get<string>("order_method", "rcm");
  typedef typename Adapter::base_adapter_t base_adapter_t;

  // TODO: Ignore case
  if (method.compare("rcm") == 0)
  {
      AlgRCM<base_adapter_t>(this->graphModel_, this->solution_, this->params_,
                      problemComm_);
  }
  else if (method.compare("Natural") == 0)
  {
      AlgNatural<base_adapter_t>(this->identifierModel_, this->solution_, this->params_, problemComm_);
  }
  else if (method.compare("Random") == 0)
  {
      AlgRandom<base_adapter_t>(this->identifierModel_, this->solution_, this->params_, problemComm_);
  }
  else if (method.compare("Minimum_Degree") == 0)
  {
      string pkg = this->params_->template get<string>("order_package", "amd");
      if (pkg.compare("amd") == 0)
      {
          AlgAMD<base_adapter_t>(this->graphModel_, this->solution_, this->params_,
                          problemComm_);
      }
  }

#ifdef HAVE_ZOLTAN2_MPI

  // The algorithm may have changed the communicator.  Change it back.

  RCP<const mpiWrapper_t > wrappedComm = rcp(new mpiWrapper_t(mpiComm_));
  problemComm_ = rcp(new Teuchos::MpiComm<int>(wrappedComm));
  problemCommConst_ = rcp_const_cast<const Comm<int> > (problemComm_);

#endif

}

////////////////////////////////////////////////////////////////////////
//template <typename Adapter>
//void OrderingProblem<Adapter>::redistribute()
//{
//  HELLO;
//}

////////////////////////////////////////////////////////////////////////
//! createOrderingProblem 
//  Method with common functionality for creating a OrderingProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void OrderingProblem<Adapter>::createOrderingProblem()
{
  HELLO;
  using Teuchos::ParameterList;

//  cout << __func__ << " input adapter type " 
//       << this->inputAdapter_->inputAdapterType() << " " 
//       << this->inputAdapter_->inputAdapterName() << endl;

#ifdef HAVE_ZOLTAN2_OVIS
  ovis_enabled(this->comm_->getRank());
#endif

  // Create a copy of the user's communicator.

  problemComm_ = this->comm_->duplicate();

#ifdef HAVE_ZOLTAN2_MPI

  // TPLs may want an MPI communicator

  Comm<int> *c = problemComm_.getRawPtr();
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm_ = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm_ = MPI_COMM_SELF;   // or would this be an error?
  }

#endif

  ParameterList *general = &(this->env_->getParametersNonConst());
  ParameterList *ordering = NULL;
  if (this->env_->hasOrderingParameters()){
    ordering = &(general->sublist("ordering"));
  }

  // Determine which parameters are relevant here.
  // For now, assume parameters similar to Zoltan:
  //   MODEL = graph, hypergraph, geometric, ids
  //   ALGORITHM = rcm, random, amd

  ModelType modelType = GraphModelType;

  typedef typename Adapter::base_adapter_t base_adapter_t;

  // Select Model based on parameters and InputAdapter type

  std::bitset<NUM_MODEL_FLAGS> graphFlags;
  std::bitset<NUM_MODEL_FLAGS> idFlags;

  switch (modelType) {

  case GraphModelType:
    graphFlags.set(SELF_EDGES_MUST_BE_REMOVED);
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, problemComm_, graphFlags));

    this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->graphModel_);

    break;



  case IdentifierModelType:
    idFlags.set(SELF_EDGES_MUST_BE_REMOVED);
    this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, problemComm_, idFlags));

    this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->identifierModel_);

    break;

  case HypergraphModelType:
  case CoordinateModelType:
    cout << __func__ << " Model type " << modelType << " not yet supported." 
         << endl;
    break;

  default:
    cout << __func__ << " Invalid model" << modelType << endl;
    break;
  }
}
} //namespace Zoltan2
#endif
