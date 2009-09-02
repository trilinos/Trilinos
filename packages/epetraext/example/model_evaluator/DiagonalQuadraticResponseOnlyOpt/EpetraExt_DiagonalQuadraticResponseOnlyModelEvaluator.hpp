#ifndef EPETRA_EXT_DIAGONAL_RESPONSE_ONLY_MODEL_EVALUATOR_HPP
#define EPETRA_EXT_DIAGONAL_RESPONSE_ONLY_MODEL_EVALUATOR_HPP


#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"


namespace EpetraExt {


/** \brief A simple quadratic parallel response-only model evaluator.
 *
 * Represents the model:
 
 \verbatim

    g[0] = 0.5 * (p-pt)^T * (p-pt)
 
 \endverbatim
 *
 * See the function <tt>evalModel()</tt> for more details.
 */
class DiagonalQuadraticResponseOnlyModelEvaluator : public EpetraExt::ModelEvaluator {
public:

  /** \brief . */
  DiagonalQuadraticResponseOnlyModelEvaluator(
    const Teuchos::RCP<Epetra_Comm> &comm,
    const int localDim, const double &pt, const double &p0, const double &scale );

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private member data

  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm_;
	Teuchos::RefCountPtr<const Epetra_Map> map_p_;
	Teuchos::RefCountPtr<const Epetra_Map> map_g_;

	double scale_;

	Teuchos::RefCountPtr<Epetra_Vector> pt_;
	Teuchos::RefCountPtr<Epetra_Vector> p0_;

  // Note defined and not to be called
  DiagonalQuadraticResponseOnlyModelEvaluator();

};


/** \brief Nonmember constructor. */
inline
Teuchos::RCP<DiagonalQuadraticResponseOnlyModelEvaluator>
diagonalQuadraticResponseOnlyModelEvaluator(
  const Teuchos::RCP<Epetra_Comm> &comm,
  const int localDim, const double &pt, const double &p0, const double &scale
  )
{
  return Teuchos::rcp(
    new DiagonalQuadraticResponseOnlyModelEvaluator(
      comm, localDim, pt, p0, scale
      )
    );
}


} // namespace EpetraExt


#endif // EPETRA_EXT_DIAGONAL_RESPONSE_ONLY_MODEL_EVALUATOR_HPP
