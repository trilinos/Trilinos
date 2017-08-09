#ifndef XPETRA_REGIONALAMG_DEF_HPP
#define XPETRA_REGIONALAMG_DEF_HPP

#include "Xpetra_RegionalAMG_decl.hpp"

namespace Xpetra{ 

	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
		RegionalAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RegionalAMG(const char* matrix_file_name, const char* elements_file_name, RCP<const Teuchos::Comm<int> > comm, Teuchos::ParameterList muelu, GlobalOrdinal num_levels, GlobalOrdinal coarsening_factor):num_levels_(num_levels),coarsening_factor_(coarsening_factor),comm_(comm), muelu_(muelu)
		{

			TEUCHOS_TEST_FOR_EXCEPT( !( matrixSplitting_.is_null() ) );
			matrixSplitting_ = rcp( new tpetra_splitting(matrix_file_name, elements_file_name, comm) );

			domainMap_ = matrixSplitting_->getDomainMap();
			rangeMap_ = matrixSplitting_->getRangeMap();
			num_regions_ = matrixSplitting_->getNumRegions();

			SetUpHierarchy();
		}

	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
		void RegionalAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetUpHierarchy()
		{
			TEUCHOS_TEST_FOR_EXCEPTION( num_levels_<=0, Exceptions::RuntimeError, "Number of levels must be a positive integer number \n" );
			TEUCHOS_TEST_FOR_EXCEPTION( num_regions_<=0 , Exceptions::RuntimeError, "No existing regions \n");

			//Creation of the MueLu list for the region multigrid preconditioner
			/*RCP<ParameterList> list = rcp(new Teuchos::ParameterList());
			list->setName("MueLu");

			list->set("verbosity", "none"); 
			list->set("number of equations", 1);
			list->set("max levels", num_levels_);
			list->set("multigrid algorithm", "unsmoothed");
			list->set("smoother: pre or post", "none");
			list->set("coarse: type", "RELAXATION"); 
			//list->set("coarse: max size", 100);

			////Geometric multigrid
			//ParameterList& factory_sublist = list->sublist("Factories");
			//ParameterList& geometric_sublist = factory_sublist.sublist("myProlongatorFact");
			//geometric_sublist.set("factory", "GeneralGeometricPFactory");
			//geometric_sublist.set("P", "myProlongatorFact");
			//geometric_sublist.set("Coarsen", "{3,3}");
			//geometric_sublist.set("order", 0);

			//Brick aggregation
			list->set("aggregation: type", "brick");
			list->set("aggregation: preserve Dirichlet nodes", true);
			list->set("aggregation: brick x size", coarsening_factor_);
			list->set("aggregation: brick y size", coarsening_factor_);
			list->set("aggregation: drop scheme", "classical");

			//Creation of Sublist for smoother	
			ParameterList& coarse_smooth_sublist = list->sublist("coarse: params");
			coarse_smooth_sublist.set("relaxation: type", "Jacobi");
			coarse_smooth_sublist.set("relaxation: sweeps", 1);
			coarse_smooth_sublist.set("relaxation: damping factor", 1.0);*/

			//ParameterList& print_sublist = list->sublist("export data");
			//print_sublist.set("P", "{1}");
			//list->print(std::cout);
			
			GlobalOrdinal num_regions = matrixSplitting_->getNumRegions();
			regionHierarchies_.resize( num_regions );

			//Instantiation of MueLu Hierarchies for each region of the domain 
			for( GlobalOrdinal region_idx = 0; region_idx<num_regions; ++region_idx )
			{

				//Creation of coordinates to pass to MueLu for the construction of the Hierarchy
				GlobalOrdinal n = matrixSplitting_->getSplittingDriver()->GetNumRegionNodes(region_idx);
				GlobalOrdinal nx = std::sqrt(n);

				RCP<const map_type> map = matrixSplitting_->getRegionMatrix(region_idx)->getRowMap();
				size_t NumMyElements = map->getNodeNumElements();
				Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

				RCP<multivector_type> coords = mv_factory_type::Build (map, 2);

				Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > Coord(2);
				Coord[0] = coords->getDataNonConst(0);
				Coord[1] = coords->getDataNonConst(1);

				//Scalar delta_x = nx / Teuchos::as<Scalar>(nx - 1);
				//Scalar delta_y = ny / Teuchos::as<Scalar>(ny - 1);

				for (size_t i = 0; i < NumMyElements; ++i) {
					GlobalOrdinal ix = MyGlobalElements[i] % nx;
					GlobalOrdinal iy = (MyGlobalElements[i] - ix) / nx;

					Coord[0][i] = Teuchos::as<Scalar>(ix);
					Coord[1][i] = Teuchos::as<Scalar>(iy);
				}

				//regionHierarchies_[region_idx] = MueLu::CreateXpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>( matrixSplitting_->getRegionMatrix( region_idx), *list, coords );
				regionHierarchies_[region_idx] = MueLu::CreateXpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>( matrixSplitting_->getRegionMatrix( region_idx), muelu_, coords );

				//Different regions may have meshes with different size, or the number of levels the user wants to create may be too big
				//with respect to the number of levels MueLu allows (i.e. the minimum coarse size may be reached for a smaller number of levels)
				//In this case, the value of num_levels_ is readjusted to the minimum mnumber of levels instantied across all the regions of the domain
				num_levels_ = std::min( num_levels_, regionHierarchies_[region_idx]->GetNumLevels() );
			}

			DefineLevels();	
	
		}


	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
		void RegionalAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DefineLevels ( )
		{ 
			Array<level> levels;
			Array<RCP<matrix_type> > P;
			P.clear();
			Array<RCP<matrix_type> > R;
			R.clear();
			Array<RCP<matrix_type> > A;
			A.clear();

			TEUCHOS_TEST_FOR_EXCEPTION( levels_.size()!=0, Exceptions::RuntimeError, "Levels structure is already initialized \n" );
			for( int i = 0; i<num_levels_; ++i )
			{
				P.clear();
				R.clear();
				A.clear();
	
				level new_level(i,num_regions_);
				for( int region_idx = 0; region_idx<num_regions_; ++region_idx )					
				{
					TEUCHOS_TEST_FOR_EXCEPTION( !regionHierarchies_[region_idx]->GetLevel(i)->IsAvailable("A") , Exceptions::RuntimeError, "No existing operator at level "<<i<<" of region "<<region_idx<<"\n");
					A.push_back( regionHierarchies_[region_idx]->GetLevel(i)->template Get<RCP<matrix_type> >("A") );
					if( i>0 )
					{
						TEUCHOS_TEST_FOR_EXCEPTION( !regionHierarchies_[region_idx]->GetLevel(i)->IsAvailable("P") , Exceptions::RuntimeError, "No existing prolongator at level "<<i<<" of region "<<region_idx<<"\n");
						TEUCHOS_TEST_FOR_EXCEPTION( !regionHierarchies_[region_idx]->GetLevel(i)->IsAvailable("R") , Exceptions::RuntimeError, "No existing restriction at level "<<i<<" of region "<<region_idx<<"\n");
						P.push_back( regionHierarchies_[region_idx]->GetLevel(i)->template Get<RCP<matrix_type> >("P") );	
						R.push_back( regionHierarchies_[region_idx]->GetLevel(i)->template Get<RCP<matrix_type> >("R") );	
					}
				}
				if( i>0 )
				{
					new_level.SetP( P );
					new_level.SetR( R );
				}
				new_level.SetA( A );

				//We pick regionToAll and we pass it to the fine level as it is
				if( 0==i )
					new_level.SetRegionToAll( matrixSplitting_->getSplittingDriver()->GetRegionToAll() );	
				else
				{
					Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > coarse_regionToAll;
					regionToAllCoarsen( levels[i-1], new_level );
				}	

				//new_level.checkConsistency();
				//new_level.ComputeRegionalJacobi();

				levels.push_back( new_level );
			}
			levels_=arcpFromArray( levels );	
		}


	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
		void RegionalAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::regionToAllCoarsen (const level& fine, level& coarse)	
		{
			TEUCHOS_TEST_FOR_EXCEPTION( fine.GetNumRegions()!=coarse.GetNumRegions(), Exceptions::RuntimeError, "Level "<<fine.GetLevelID()<<"has "<<fine.GetNumRegions()<<" regions instantiated whereas level "<<coarse.GetLevelID()<<"has "<<coarse.GetNumRegions()<<" regions instantiated \n" );

			Array<Array<std::tuple<GlobalOrdinal,GlobalOrdinal> > > fine_regionToAll = fine.GetRegionToAll();
			Array<Array<std::tuple<GlobalOrdinal,GlobalOrdinal> > > coarse_regionToAll(num_regions_);

			for( GlobalOrdinal region_idx = 0; region_idx<num_regions_; ++region_idx )
			{

				//Coordinates needed to apply the coarsening
				GlobalOrdinal n = fine.GetNumRegionNodes(region_idx);
				GlobalOrdinal nx = std::sqrt(n);
				GlobalOrdinal ny = nx;

				TEUCHOS_TEST_FOR_EXCEPTION( (nx-1)%coarsening_factor_!=0 || (ny-1)%coarsening_factor_!=0, Exceptions::RuntimeError, "Region: "<<region_idx<<" cannot be coarsened by factor equal to "<<coarsening_factor_<<"n: "<<n<<" - nx: "<<nx<<" - ny: "<<ny<<" \n" );

				//here below we explot the fact that the regionToAll has been sorted in ascending order for the region node index
				GlobalOrdinal count = 1;
				
				while( count<=fine_regionToAll[region_idx].size() )
				{
					coarse_regionToAll[region_idx].push_back(fine_regionToAll[region_idx][count]);
					if( count%ny!=0 )
						count = count + coarsening_factor_;
					else
						count = count + ( coarsening_factor_-1 ) * ny + 1;

				}

				coarse.SetRegionToAll( coarse_regionToAll );

			}

		}



	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
		void
			RegionalAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply (const multivector_type& X, multivector_type& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta)const
		{ 


		}

}
#endif
