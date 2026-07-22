// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_ColoringProblem.hpp
    \brief Defines the ColoringProblem class.
*/

#ifndef _ZOLTAN2_COLORINGPROBLEM_HPP_
#define _ZOLTAN2_COLORINGPROBLEM_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_ColoringAlgorithms.hpp>
#include <Zoltan2_ColoringSolution.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <string>

#include <bitset>

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////

/*! \brief ColoringProblem sets up coloring problems for the user.
 *
 *  The ColoringProblem is the core of the Zoltan2 coloring API.
 *  Based on the the user's input and parameters, the ColoringProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the ColoringProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo - Should Problems and Solution have interfaces for returning
 *          views and for returning RCPs?  Or just one?  At a minimum,
 *          we should have the word "View" in function names that return views.
 *
 *  \todo - Currently, only serial and shared-memory coloring are supported.
 */

template<typename Adapter>
class ColoringProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

#ifdef HAVE_ZOLTAN2_MPI
   typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
#endif

  /*! \brief Destructor
   */
  virtual ~ColoringProblem() {};

  /*! \brief Constructor that uses a Teuchos::Comm
   */
  ColoringProblem(Adapter *A, ParameterList *p,
                  const Teuchos::RCP<const Teuchos::Comm<int> > &comm) :
    Problem<Adapter>(A, p, comm)
  {
    HELLO;
    createColoringProblem();
  };

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor that takes an MPI communicator
   */
  ColoringProblem(Adapter *A, ParameterList *p, MPI_Comm mpicomm) :
  ColoringProblem(A, p,
                  rcp<const Comm<int> >(new Teuchos::MpiComm<int>(
                                            Teuchos::opaqueWrapper(mpicomm))))
  {}
#endif

  /*! \brief Constructor that uses a default communicator
   */
  ColoringProblem(Adapter *A, ParameterList *p) :
  ColoringProblem(A, p, Tpetra::getDefaultComm())
  {}

  /*! \brief Set up validators specific to this Problem
  */
  static void getValidParameters(ParameterList & pl)
  {
    RCP<Teuchos::StringValidator> color_method_Validator = Teuchos::rcp(
      new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "SerialGreedy","D1","D1-2GL","D2","PD2" )));
    pl.set("color_method", "SerialGreedy", "coloring algorithm",
     color_method_Validator);
    pl.set("verbose", false, "print all output", Environment::getBoolValidator());
    pl.set("timing", false, "print timing data", Environment::getBoolValidator());
    pl.set("serial_threshold",0,"vertices to recolor in serial",Environment::getAnyIntValidator());
    pl.set("recolor_degrees",true,"recolor based on vertex degrees",Environment::getBoolValidator());
  }

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

  ColoringSolution<Adapter> *getSolution() {
    // Get the raw ptr from the rcp
    return solution_.getRawPtr();
  };

private:
  void createColoringProblem();

  RCP<ColoringSolution<Adapter> > solution_;
  size_t localNumObjects_;
};


////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void ColoringProblem<Adapter>::solve(bool newData)
{
  HELLO;

  try
  {
      this->solution_ = rcp(new ColoringSolution<Adapter>(localNumObjects_));
  }
  Z2_FORWARD_EXCEPTIONS;

  // Determine which algorithm to use based on defaults and parameters.
  // Need some exception handling here, too.

  std::string method = this->params_->template get<std::string>("color_method", "SerialGreedy");

  try
  {
  // TODO: Ignore case
  if (method.compare("SerialGreedy") == 0)
  {
      modelFlag_t graphFlags;
      graphFlags.set(REMOVE_SELF_EDGES);
      graphFlags.set(BUILD_LOCAL_GRAPH);

      AlgSerialGreedy<Adapter> alg(this->inputAdapter_, this->params_,
                                   this->env_, this->comm_, graphFlags);
      alg.color(this->solution_);
  }
  else if (method.compare("D1") == 0)
  {
      AlgDistance1<Adapter> alg(this->inputAdapter_, this->params_,
		                this->env_, this->comm_);
      alg.color(this->solution_);
  }
  else if (method.compare("D1-2GL") == 0)
  {
      AlgDistance1TwoGhostLayer<Adapter> alg(this->inputAdapter_,this->params_,
		                             this->env_, this->comm_);
      alg.color(this->solution_);
  } else if(method.compare("D2") == 0)
  {
      AlgDistance2<Adapter> alg(this->inputAdapter_, this->params_,
 		                this->env_, this->comm_);
      alg.color(this->solution_);
  } else if (method.compare("PD2") == 0)
  {
      AlgPartialDistance2<Adapter> alg(this->inputAdapter_, this->params_,
		                       this->env_, this->comm_);
      alg.color(this->solution_);
  }
  }
  Z2_FORWARD_EXCEPTIONS;
}

////////////////////////////////////////////////////////////////////////
//template <typename Adapter>
//void ColoringProblem<Adapter>::redistribute()
//{
//  HELLO;
//}

////////////////////////////////////////////////////////////////////////
//! createColoringProblem
//  Method with common functionality for creating a ColoringProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void ColoringProblem<Adapter>::createColoringProblem()
{
  HELLO;
  using Teuchos::ParameterList;

//  std::cout << __func__zoltan2__ << " input adapter type "
//       << this->inputAdapter_->inputAdapterType() << " "
//       << this->inputAdapter_->inputAdapterName() << std::endl;

  // Create a copy of the user's communicator.

  // Only graph model supported.
  // TODO: Allow hypergraph later?

  ModelType modelType = GraphModelType;
  const auto adapterType = this->baseInputAdapter_->adapterType();

  // Select Model based on parameters and InputAdapter type

  switch (modelType)
  {

  case GraphModelType:
  {
    switch (adapterType)
    {
    case MatrixAdapterType:
    {
      localNumObjects_ = this->baseInputAdapter_->getLocalNumIDs();
    }
    break;

    case GraphAdapterType:
    {
      const auto ia = dynamic_cast<const GraphAdapter<user_t> *>(&(*(this->baseInputAdapter_)));
      localNumObjects_ = ia->getLocalNumVertices();
    }
    break;

    case MeshAdapterType:
    {
      const auto ia = dynamic_cast<const MeshAdapter<user_t> *>(&(*(this->baseInputAdapter_)));
      localNumObjects_ = ia->getLocalNumOf(ia->getPrimaryEntityType());
    }
    break;

    default:
    {
      // Avoid warning
    }
    }
  }
  break;

  case IdentifierModelType:
  case HypergraphModelType:
  case CoordinateModelType:
    std::cout << __func__zoltan2__ << " Model type " << modelType
              << " not yet supported." << std::endl;
    break;

  default:
    std::cout << __func__zoltan2__ << " Invalid model" << modelType
              << std::endl;
    break;
  }
}
} //namespace Zoltan2

#endif
