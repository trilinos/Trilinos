// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_TPETRA_ME_1DFEM_DECL_HPP
#define NOX_TPETRA_ME_1DFEM_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_TimeMonitor.hpp"

template<class Scalar, class LO, class GO, class Node>
class EvaluatorTpetra1DFEM;

/** \brief Nonmember constuctor.
 *
 * \relates EvaluatorTpetra1DFEM
 */
template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar, LO, GO, Node> >
evaluatorTpetra1DFEM(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Tpetra::global_size_t numGlobalElements,
                     const Scalar zMin,
                     const Scalar zMax);


/** \brief 1D Finite Element model for nonlinear heat conduction
 *
 * The equation modeled is:

 \verbatim

   d2T
   --- - k * T**2 = 0
   dz2

   subject to:
      T  = p(4), defaults to 1.0 @ z = zMin
      T' = 0.0 @ z = zMax
      k = p(2), defaults to 1.0 (independent parameter)

   For LOCA testing, parameters and responses are used to evaluate a
   contraint equation. The parameter k is the unknown variable used to
   enforce the response. We add in dummy responses for unit testing
   correct offsets in LOCA parameter handling code.:

      p(0) is a dummy parameter (throws if queried)
      p(1) is a dummy parameter (throws if queried)
      p(2) is a constant "k" multiplier on the source term
      p(3) is a dummy parameter (throws if queried)
      p(4) is the Dirichlet BC temperature value for left side node "T_left"

      g(0) is a dummy response (throws if queried)
      g(1) is a dummy response (throws if queried)
      g(2) is a dummy response (throws if queried)
      g(3) is a dummy response (throws if queried)
      g(4) = T(zMax) - 2.0
      g(5) is a dummy response (throws if queried)
      g(6) = 2.0 * T(zMin) - T(zMax)

 \endverbatim

 * The Matrix <tt>W = d(f)/d(x)</tt> is implemented as a
 * <tt>Thyra::LinearOpBase</tt> object and the class
 * <tt>Thyra::DefaultSerialDenseLinearOpWithSolveFactory</tt> is used to
 * create the linear solver.
 */
template<class Scalar, class LO, class GO, class Node>
class EvaluatorTpetra1DFEM : public ::Thyra::ModelEvaluator<Scalar>
{
public:

  // Public typedefs
  typedef Scalar scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef Node node_type;
  typedef Tpetra::Map<LO, GO, Node> tpetra_map;
  typedef Tpetra::CrsGraph<LO, GO, Node> tpetra_graph;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> tpetra_matrix;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> tpetra_vec;
  typedef ::Thyra::VectorBase<Scalar> thyra_vec;
  typedef ::Thyra::MultiVectorBase<Scalar> thyra_mvec;
  typedef ::Thyra::VectorSpaceBase<Scalar> thyra_vec_space;
  typedef ::Thyra::LinearOpBase<Scalar> thyra_op;
  typedef ::Thyra::PreconditionerBase<Scalar> thyra_prec;

  // Constructor
  EvaluatorTpetra1DFEM(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                       const Tpetra::global_size_t numGlobalElements,
                       const Scalar zMin,
                       const Scalar zMax);

  /** \name Initializers/Accessors */
  //@{

  /** \brief . */
  void set_x0(const Teuchos::ArrayView<const Scalar> &x0);

  /** \brief . */
  void setShowGetInvalidArgs(bool showGetInvalidArg);

  void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  Teuchos::RCP<const thyra_vec_space> get_x_space() const;
  Teuchos::RCP<const thyra_vec_space> get_f_space() const;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<thyra_op> create_W_op() const;
  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  Teuchos::RCP<thyra_prec> create_W_prec() const;

  // These are for constraint solver support
  int Np () const;
  int Ng () const;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  Teuchos::ArrayView<const std::string> get_g_names(int j) const;
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_DfDp_op(int l) const;
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar> > create_DgDx_op(int j) const;
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar> > create_DgDx_dot_op(int j) const;
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const;

  void evalModel(const ::Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                 const ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  // Dummy
  Teuchos::RCP<const thyra_vec_space> get_f_multiplier_space() const{return Teuchos::null;}
  Teuchos::RCP<const thyra_vec_space> get_g_multiplier_space(int j) const{return Teuchos::null;}
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const{return prototypeInArgs_;}
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const{return prototypeInArgs_;}
  Teuchos::RCP<::Thyra::LinearOpWithSolveBase<Scalar>> create_W() const{return Teuchos::null;}
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_hess_f_xx() const{return Teuchos::null;}
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_hess_f_xp(int l) const{return Teuchos::null;}
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_hess_f_pp( int l1, int l2 ) const{return Teuchos::null;}
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_hess_g_xx(int j) const{return Teuchos::null;}
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_hess_g_xp( int j, int l ) const{return Teuchos::null;}
  Teuchos::RCP<::Thyra::LinearOpBase<Scalar>> create_hess_g_pp( int j, int l1, int l2 ) const{return Teuchos::null;}
  void reportFinalPoint (const ::Thyra::ModelEvaluatorBase::InArgs<Scalar> &finalPoint, const bool wasSolved);
  //@}

  std::pair<int,int> get_p_index(const std::string& p_name) const;
  std::pair<int,int> get_g_index(const std::string& g_name) const;

private:

  /** Allocates and returns the Jacobian matrix graph */
  virtual Teuchos::RCP<const tpetra_graph> createGraph();

private: // data members

  const Teuchos::RCP<const Teuchos::Comm<int> >  comm_;
  const Tpetra::global_size_t numGlobalElements_;
  const Scalar zMin_;
  const Scalar zMax_;
  const int Np_; // Number of parameters
  const int Ng_; // Number of responses
  const bool printDebug_;

  Teuchos::RCP<const thyra_vec_space> xSpace_;
  Teuchos::RCP<const tpetra_map>   xOwnedMap_;
  Teuchos::RCP<const tpetra_map>   xGhostedMap_;
  Teuchos::RCP<const Tpetra::Import<LO, GO, Node> > importer_;

  Teuchos::RCP<const thyra_vec_space> fSpace_;
  Teuchos::RCP<const tpetra_map>   fOwnedMap_;

  Teuchos::RCP<const tpetra_graph>  W_graph_;

  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

  Teuchos::RCP<tpetra_vec> nodeCoordinates_;

  mutable Teuchos::RCP<tpetra_vec> uPtr_;
  mutable Teuchos::RCP<tpetra_vec> xPtr_;

  mutable Teuchos::RCP<tpetra_vec> J_diagonal_;

  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<thyra_vec> x0_;
  bool showGetInvalidArg_;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  mutable Teuchos::RCP<Teuchos::Time> residTimer_;
  mutable Teuchos::RCP<Teuchos::Time> jacTimer_;

  // Optional parameter and response support for LOCA
  std::vector<Teuchos::RCP<Teuchos::Array<std::string>>> pNames_;
  std::vector<std::vector<std::string>> gNames_;
  Teuchos::RCP<const tpetra_map> pMap_; // locally replicated scalar
  Teuchos::RCP<const thyra_vec_space> pSpace_;
  Teuchos::RCP<thyra_vec> p2_; // k value in equation
  Teuchos::RCP<thyra_vec> p4_; // T_left value in equation Dirichlet BC
  Teuchos::RCP<const tpetra_map> gMap_; // locally replicated scalar
  Teuchos::RCP<const thyra_vec_space> gSpace_;
  Teuchos::RCP<const tpetra_map> dgdpMap_; // locally replicated dense matrix
  Teuchos::RCP<const thyra_vec_space> dgdpSpace_;

  std::unordered_map<std::string,std::pair<int,int>> p_name_to_index_;
  std::unordered_map<std::string,std::pair<int,int>> g_name_to_index_;
};

//==================================================================
#include "ME_Tpetra_1DFEM_def.hpp"
//==================================================================

#endif
