#ifndef PANZER_MODEL_EVALUATOR_HPP
#define PANZER_MODEL_EVALUATOR_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"

#include <Kokkos_DefaultNode.hpp>

namespace panzer {

  template<typename Scalar, typename LO, typename GO, typename NODE>
class ModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

//   typedef typename panzer::Traits<T>::lid_type LO;
//   typedef typename panzer::Traits<T>::gid_type GO;
//   typedef typename panzer::Traits<T>::node_type NODE;

public:

  /** \name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  ModelEvaluator();

  /** \brief . */
  void set_d(const Scalar &d);

  /** \brief . */
  void set_p(const Teuchos::ArrayView<const Scalar> &p);

  /** \brief . */
  void set_x0(const Teuchos::ArrayView<const Scalar> &x0);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private: // data members

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Scalar d_;
  Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO,NODE> > x0_;
  Teuchos::Array<Scalar> p_;
  Teuchos::RCP<const Tpetra::CrsGraph<LO,GO,NODE> > W_op_graph_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;
};


}

#include "Panzer_ModelEvaluatorT.hpp"

#endif 
