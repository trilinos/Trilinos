#ifndef XPETRA_REGIONALAMG_DECL_HPP
#define XPETRA_REGIONALAMG_DECL_HPP

//Xpetra
#include <Xpetra_Operator.hpp>
#include "Xpetra_MatrixSplitting.hpp"
#include "Xpetra_Level_def.hpp"

//MueLu
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>


//This class inherits from Xpetra::Operator because the goal is to allow its usage 
//as a preconditioner for a Belos linear solver

namespace Xpetra{

	template<class Scalar					= MultiVector<>::scalar_type,
		class LocalOrdinal   				= typename MultiVector<Scalar>::local_ordinal_type,
		class GlobalOrdinal  				= typename MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
		class Node           				= typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
		class RegionalAMG : Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {

		typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
		typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal> multivector_type;
		typedef Matrix< Scalar, LocalOrdinal, GlobalOrdinal, Node > matrix_type;
		typedef MatrixSplitting<Scalar,LocalOrdinal,GlobalOrdinal,Node,Xpetra::UseTpetra, false> tpetra_splitting;
		typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> Hierarchy;
		typedef MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> mv_factory_type;
		typedef Level<Scalar,LocalOrdinal,GlobalOrdinal,Node> level;

		public: 

			//! Constructors
			//@{

			RegionalAMG()
			{
				std::cout<<"This version of constructor is not implemented yet"<<std::endl;
			};
			
			//! This is the actual constructor we are interested in calling
			RegionalAMG( const char* , const char* , RCP<const Teuchos::Comm<int> >, Teuchos::ParameterList, GlobalOrdinal, GlobalOrdinal );

			//@}

			GlobalOrdinal GetNumLevels(){return num_levels_;}

			//! Methods to extract Map information
			//@{
			//! For now, the domain Map coincides with the Domain Map of the composite matrix at the fine level
			virtual RCP<const map_type > getDomainMap()const{return domainMap_;} 

			//! For now, the domain Map coincides with the Range Map of the composite matrix at the fine level
			virtual RCP<const map_type > getRangeMap()const{return rangeMap_;}

			//@}
			
			//! Apply method
			//@{
			//!N.B.: Still to implement
			virtual void
				apply (const multivector_type& X, multivector_type& Y, 
				Teuchos::ETransp mode = Teuchos::NO_TRANS, 
				Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(), 
				Scalar beta = Teuchos::ScalarTraits<Scalar>::zero())const; 
			//@}

			//! hasTransposeApply should not be needed
			virtual bool hasTransposeApply() const { return false; }

		private:

			//Total number of levels in the hierarchy
			GlobalOrdinal num_levels_ = -1;

			//Total number of regions
			GlobalOrdinal num_regions_ = -1;

			//Coarsening factor to transfer quantities across levels
			int coarsening_factor_ = -1;

			RCP<const Teuchos::Comm<int> > comm_;

			Teuchos::ParameterList muelu_;

			RCP<const map_type > domainMap_;
			RCP<const map_type > rangeMap_;

			//matrixSplitting associated with the composite matrix at the fine level
			RCP<tpetra_splitting> matrixSplitting_;

			//Array of MueLu hierarchies (one for each region)
			Array<RCP<Hierarchy> > regionHierarchies_;

			//Array of levels in the new hierarchy (each levels contains quantities associaed with every region for that level)
			Array<RCP<level> > levels_;

			//This methods construct the levels using the quantities stored in the MueLu hierarchies
			void SetUpHierarchy();

			void DefineLevels();

			virtual void regionToAllCoarsen(const level&, level& );

	};

} //namespace Xpetra

#endif
