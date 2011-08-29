#ifndef PIRO_EPETRA_NE_COUPLED_MODEL_EVALUATOR_HPP
#define PIRO_EPETRA_NE_COUPLED_MODEL_EVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Piro_Epetra_StokhosSolver.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LocalMap.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "EpetraExt_MultiComm.h"

namespace Piro {
  namespace Epetra {

    class NECoupledModelEvaluator : public EpetraExt::ModelEvaluator {
    public:

      /** \brief . */
      NECoupledModelEvaluator(
	const Teuchos::RCP<Teuchos::ParameterList>& piroParams,
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelA, 
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelB,
	const Teuchos::RCP<Teuchos::ParameterList>& piroParamsA,
	const Teuchos::RCP<Teuchos::ParameterList>& piroParamsB,
	const Teuchos::RCP<const Epetra_Comm>& comm);

      /** \name Overridden from EpetraExt::ModelEvaluator . */
      //@{

      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
      //! Return array of parameter names
      Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
      /** \brief . */
      Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
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
      
      typedef Stokhos::StandardStorage<int,double> StorageType;

      Teuchos::RCP<EpetraExt::ModelEvaluator> modelA;
      Teuchos::RCP<EpetraExt::ModelEvaluator> modelB;
      Teuchos::RCP<Teuchos::ParameterList> piroParamsA;
      Teuchos::RCP<Teuchos::ParameterList> piroParamsB;
      Teuchos::RCP<Teuchos::ParameterList> piroParams;
      Teuchos::RCP<const Epetra_Comm> comm;

      Teuchos::RCP<EpetraExt::ModelEvaluator> solverA;
      Teuchos::RCP<EpetraExt::ModelEvaluator> solverB;
      Teuchos::RCP<Piro::Epetra::StokhosSolver> sgSolverA;
      Teuchos::RCP<Piro::Epetra::StokhosSolver> sgSolverB;

      Teuchos::RCP<const Epetra_LocalMap> map_p;
      Teuchos::RCP<const Epetra_LocalMap> map_rvar;
      Teuchos::RCP<const Epetra_LocalMap> map_g;
      Teuchos::RCP<const Epetra_LocalMap> map_x;
      Teuchos::RCP<const Epetra_Map> x_OverlapMap;
      Teuchos::RCP<const Epetra_Map> g_OverlapMap;
      Teuchos::RCP<Epetra_Vector> x0;
      Teuchos::RCP<Epetra_Vector> p1;
      Teuchos::RCP<Epetra_Vector> p2;
      Teuchos::RCP<Epetra_Vector> rVars;
      Teuchos::RCP<Epetra_Import> OtoLx;
      Teuchos::RCP<Epetra_Import> LtoOx;
      Teuchos::RCP<Epetra_Import> OtoLg;
      Teuchos::RCP<Epetra_Import> LtoOg;
      Teuchos::RCP<Epetra_MultiVector> dirvec1;
      Teuchos::RCP<Epetra_MultiVector> dirvec2;
      Teuchos::RCP<Epetra_Vector> g1;
      Teuchos::RCP<Epetra_Vector> g2;
      
      Teuchos::RCP<EpetraExt::ModelEvaluator::Derivative> dgdp1;
      Teuchos::RCP<EpetraExt::ModelEvaluator::Derivative> dgdp2;
      Teuchos::RCP<Teuchos::Array<std::string> > pnames;
      mutable Teuchos::RCP<Epetra_CrsMatrix> A;
      
      mutable Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis;
      mutable Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
      mutable Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad;
      mutable Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,StorageType> > expansion;
      mutable Teuchos::RCP<const Epetra_BlockMap> sg_overlap_map;
      mutable OutArgs::sg_vector_t p1_sg;
      mutable OutArgs::sg_vector_t p2_sg;
      mutable OutArgs::sg_vector_t rVars_sg;
      mutable OutArgs::sg_vector_t g1_sg;
      mutable OutArgs::sg_vector_t g2_sg;
      mutable Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp1_sg;
      mutable Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp2_sg;
      
      mutable Teuchos::RCP<EpetraExt::ModelEvaluator> DiffProbLocal;
      mutable Teuchos::RCP<EpetraExt::ModelEvaluator> HeatProbLocal;

      int MyPID;
      int NumProc;
      int DnumRandVar;
      int HnumRandVar;
      int totNumRandVar;
      bool reduce_dimension;
      bool orthogonalize_bases;
      bool eval_W_with_f;
    };

  }

}

#endif
