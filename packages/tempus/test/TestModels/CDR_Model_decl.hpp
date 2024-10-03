//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_CDR_MODEL_DECL_HPP
#define TEMPUS_CDR_MODEL_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

class Epetra_Comm;
class Epetra_Map;
class Epetra_Vector;
class Epetra_CrsGraph;
class Epetra_Import;

namespace Tempus_Test {

template <class Scalar>
class ModelEvaluator1DFEM;

/** \brief 1D CGFEM model for convection/diffusion/reaction
 *
 * The equation modeled is:

 \verbatim

dT     dT   d^2(T)
-- + a -- + ------ - K * T**2 = 0
dt     dz    dz^2

   subject to:
      T = 1.0 @ z = z_min

 \endverbatim

 * The Matrix <tt>W = d(f)/d(x)</tt> is implemented as a
 * <tt>Thyra::MultiVectorBase</tt> object and the class
 * <tt>Thyra::DefaultSerialDenseLinearOpWithSolveFactory</tt> is used to
 * create the linear solver.
 */
template <class Scalar>
class CDR_Model : public ::Thyra::StateFuncModelEvaluatorBase<Scalar> {
 public:
  CDR_Model(const Teuchos::RCP<const Epetra_Comm> &comm,
            const int num_global_elements, const Scalar z_min,
            const Scalar z_max,
            const Scalar a,   // convection
            const Scalar k);  // source

  /** \name Initializers/Accessors */
  //@{

  void set_x0(const Teuchos::ArrayView<const Scalar> &x0);

  void setShowGetInvalidArgs(bool showGetInvalidArg);

  void set_W_factory(
      const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >
          &W_factory);

  //@}

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > create_W() const;
  Teuchos::RCP< ::Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar> > create_W_prec() const;
  //@}

 private:
  /** Allocates and returns the Jacobian matrix graph */
  virtual Teuchos::RCP<Epetra_CrsGraph> createGraph();

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{

  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
      const ::Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
      const ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  //@}

 private:  // data members
  const Teuchos::RCP<const Epetra_Comm> comm_;
  const int num_global_elements_;
  const Scalar z_min_;
  const Scalar z_max_;
  const Scalar a_;
  const Scalar k_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Epetra_Map> x_owned_map_;
  Teuchos::RCP<const Epetra_Map> x_ghosted_map_;
  Teuchos::RCP<const Epetra_Import> importer_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Epetra_Map> f_owned_map_;

  Teuchos::RCP<Epetra_CrsGraph> W_graph_;

  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

  Teuchos::RCP<Epetra_Vector> node_coordinates_;
  Teuchos::RCP<Epetra_Vector> ghosted_node_coordinates_;

  mutable Teuchos::RCP<Epetra_Vector> u_ptr;
  mutable Teuchos::RCP<Epetra_Vector> u_dot_ptr;
  mutable Teuchos::RCP<Epetra_Vector> x_ptr;

  mutable Teuchos::RCP<Epetra_Vector> J_diagonal_;

  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > x0_;
  Teuchos::Array<Scalar> p_;
  bool showGetInvalidArg_;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;
};

//==================================================================
// Finite Element Basis Object
class Basis {
 public:
  // Constructor
  Basis();

  // Destructor
  ~Basis();

  // Calculates the values of u and x at the specified gauss point
  void computeBasis(int gp, double *z, double *u, double *u_dot = nullptr);

 public:
  // Variables that are calculated at the gauss point
  double *phi, *dphide;
  double uu, zz, duu, eta, wt;
  double dz;
  // These are only needed for transient
  double uu_dot, duu_dot;
};

}  // namespace Tempus_Test

#endif
