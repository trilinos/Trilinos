#include "Smooth_Prolongation.hpp"

namespace MueLu {

AdditiveVariant::AdditiveVariant(RCP<crs_matrix_type> A, RCP<multivector_type> coords, DomainPartitioning domain) {
  domain_ = domain;

  TEUCHOS_TEST_FOR_EXCEPT(!(GlobalComm_.is_null()));
  GlobalComm_ = A->getComm();
  TEUCHOS_TEST_FOR_EXCEPT(GlobalComm_.is_null());

  TEUCHOS_TEST_FOR_EXCEPT(!(coords_.is_null()));
  coords_ = coords;
  TEUCHOS_TEST_FOR_EXCEPT((coords_.is_null()));

  TEUCHOS_TEST_FOR_EXCEPT(!(DomainMap_.is_null()));
  DomainMap_ = A->getDomainMap();
  TEUCHOS_TEST_FOR_EXCEPT(DomainMap_.is_null());

  TEUCHOS_TEST_FOR_EXCEPT(!(RangeMap_.is_null()));
  RangeMap_ = A->getRangeMap();
  TEUCHOS_TEST_FOR_EXCEPT(RangeMap_.is_null());

  TEUCHOS_TEST_FOR_EXCEPT(!(B_fine_.is_null()) || !(B_coarse_.is_null()));
  AdditiveFineSmoother(A);
  AdditiveCoarseSolver(A);
  TEUCHOS_TEST_FOR_EXCEPT(B_fine_.is_null() || B_coarse_.is_null());

  /*TEUCHOS_TEST_FOR_EXCEPT( !( B_DD_.is_null() ) || !( B_coarse_.is_null() ) );
  AdditiveFineSmoother( A );
  AdditiveCoarseSolver( A );
  TEUCHOS_TEST_FOR_EXCEPT( B_DD_.is_null() || B_coarse_.is_null() );
  TEUCHOS_TEST_FOR_EXCEPT( !B_DD_->isInitialized() || !B_DD_->isComputed() );*/
}

void AdditiveVariant::AdditiveFineSmoother(RCP<crs_matrix_type> A) {
  // Creation of the MueLu list for the DD preconditioner
  RCP<ParameterList> dd_list = rcp(new Teuchos::ParameterList());
  dd_list->setName("MueLu");
  dd_list->set("verbosity", "low");
  dd_list->set("number of equations", 1);
  dd_list->set("max levels", 1);
  dd_list->set("coarse: type", "SCHWARZ");  // FOR A ONE LEVEL PRECONDITIONER THE COARSE LEVEL IS INTERPRETED AS SMOOTHING LEVEL

  ParameterList& dd_smooth_sublist = dd_list->sublist("coarse: params");
  dd_smooth_sublist.set("schwarz: overlap level", 0);
  dd_smooth_sublist.set("schwarz: combine mode", "Zero");
  dd_smooth_sublist.set("subdomain solver name", "RILUK");

  ParameterList& coarse_subdomain_solver = dd_smooth_sublist.sublist("subdomain solver parameters");
  coarse_subdomain_solver.set("fact: iluk level-of-fill", 3);
  coarse_subdomain_solver.set("fact: absolute threshold", 0.);
  coarse_subdomain_solver.set("fact: relative threshold", 1.);
  coarse_subdomain_solver.set("fact: relax value", 0.);

  B_fine_ = CreateTpetraPreconditioner((RCP<operator_type>)A, *dd_list);

  /*B_DD_ = rcp( new Ifpack2::AdditiveSchwarz< row_matrix_type, precond_type >(A, 0) );
  B_DD_->initialize();
  B_DD_->compute();*/
}

void AdditiveVariant::AdditiveCoarseSolver(RCP<crs_matrix_type> A) {
  // Creation of the MueLu list for the DD preconditioner
  RCP<ParameterList> coarse_list = rcp(new Teuchos::ParameterList());
  coarse_list->setName("MueLu");
  coarse_list->set("verbosity", "low");
  coarse_list->set("number of equations", 1);
  coarse_list->set("max levels", 2);
  coarse_list->set("multigrid algorithm", "unsmoothed");
  coarse_list->set("aggregation: type", "brick");
  coarse_list->set("aggregation: brick x size", domain_.bricksize_x);
  coarse_list->set("aggregation: brick y size", domain_.bricksize_y);
  coarse_list->set("aggregation: brick z size", domain_.bricksize_z);
  coarse_list->set("aggregation: drop scheme", "classical");
  coarse_list->set("smoother: pre or post", "none");
  coarse_list->set("repartition: enable", true);
  coarse_list->set("repartition: partitioner", "zoltan");
  coarse_list->set("repartition: start level", 1);
  coarse_list->set("repartition: min rows per proc", static_cast<int>(A->getGlobalNumRows()));
  coarse_list->set("repartition: max imbalance", 1.2);
  coarse_list->set("repartition: remap parts", false);
  coarse_list->set("coarse: type", "SCHWARZ");

  // Creation of Sublist for smoother
  ParameterList& coarse_smooth_sublist = coarse_list->sublist("coarse: params");
  coarse_smooth_sublist.set("schwarz: overlap level", 0);
  coarse_smooth_sublist.set("schwarz: combine mode", "Zero");
  coarse_smooth_sublist.set("subdomain solver name", "RILUK");

  // Creation of the sublist for the subdomain solver
  ParameterList& coarse_subdomain_solver = coarse_smooth_sublist.sublist("subdomain solver parameters");
  coarse_subdomain_solver.set("fact: iluk level-of-fill", 3);
  coarse_subdomain_solver.set("fact: absolute threshold", 0.);
  coarse_subdomain_solver.set("fact: relative threshold", 1.);
  coarse_subdomain_solver.set("fact: relax value", 0.);

  /*ParameterList& coarse_file_sublist = coarse_list->sublist("export data");
coarse_file_sublist.set("A", "{0}");
coarse_file_sublist.set("P", "{0,1}");
coarse_file_sublist.set("R", "{1}");*/

  // Manual set up of the prolongation and restriction
  MueLu::ParameterListInterpreter<scalar_type> mueLuFactory(*coarse_list);
  RCP<MueLu::Hierarchy<scalar_type>> H = mueLuFactory.CreateHierarchy();
  H->setVerbLevel(Teuchos::VERB_HIGH);
  RCP<xpetra_matrix> mueluA = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(A);

  RCP<Xpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>> coords = Xpetra::toXpetra(coords_);

  H->GetLevel(0)->Set("A", mueluA);
  H->GetLevel(0)->Set("Coordinates", coords);

  // Multigrid setup phase
  mueLuFactory.SetupHierarchy(*H);

  RCP<Level> L = H->GetLevel(1);

  RCP<Xpetra::Matrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>> prolong, restr;

  if (L->IsAvailable("P"))
    prolong = L->template Get<RCP<Xpetra::Matrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>>>("P");

  if (L->IsAvailable("R"))
    restr = L->template Get<RCP<Xpetra::Matrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>>>("R");

  RCP<crs_matrix_type> tpetra_prolong = MueLuUtilities::Op2NonConstTpetraCrs(prolong);
  RCP<crs_matrix_type> tpetra_restr   = MueLuUtilities::Op2NonConstTpetraCrs(restr);

  int mypid = GlobalComm_->getRank();
  GlobalComm_->barrier();

  // We have to transform P into a condensed multivector
  RCP<multivector_type> identity_shrunk                                  = rcp(new multivector_type(tpetra_prolong->getDomainMap(), 3));
  Teuchos::ArrayView<const global_ordinal_type> myIdentityGlobalElements = tpetra_prolong->getDomainMap()->getLocalElementList();
  typedef typename Teuchos::ArrayView<const global_ordinal_type>::const_iterator iter_type;

  int my_color                           = (mypid - 1) % 3;
  Teuchos::ArrayRCP<scalar_type> localMV = identity_shrunk->getDataNonConst(my_color);

  for (iter_type it = myIdentityGlobalElements.begin(); it != myIdentityGlobalElements.end(); ++it) {
    const local_ordinal_type i_local = *it;
    const local_ordinal_type aux     = identity_shrunk->getMap()->getLocalElement(i_local);
    localMV[aux]                     = 1.0;
  }

  RCP<multivector_type> P_shrunk = rcp(new multivector_type(tpetra_prolong->getRangeMap(), 3));
  tpetra_prolong->apply(*identity_shrunk, *P_shrunk);
  RCP<multivector_type> AP_shrunk = rcp(new multivector_type(A->getRangeMap(), 3));
  A->apply(*P_shrunk, *AP_shrunk);

  TEUCHOS_TEST_FOR_EXCEPT(B_fine_.is_null());

  //========================================================================================================

  // CREATION OF BAP

  RCP<multivector_type> BAP_multivector = rcp(new multivector_type(B_fine_->getRangeMap(), AP_shrunk->getNumVectors()));
  RCP<multivector_type> BAP_shrunk      = rcp(new multivector_type(B_fine_->getRangeMap(), AP_shrunk->getNumVectors()));
  B_fine_->apply(*AP_shrunk, *BAP_shrunk, Teuchos::NO_TRANS, Teuchos::ScalarTraits<scalar_type>::one(), Teuchos::ScalarTraits<scalar_type>::zero());

  // I just need this to generate the right colMap to populate BAP
  RCP<crs_matrix_type> AP = rcp(new crs_matrix_type(tpetra_prolong->getRowMap(), tpetra_prolong->getColMap(), tpetra_prolong->getGlobalNumRows()));
  Tpetra::MatrixMatrix::Multiply(*A, false, *tpetra_prolong, false, *AP, true);

  GlobalComm_->barrier();

  RCP<crs_matrix_type> BAP                                      = rcp(new crs_matrix_type(tpetra_prolong->getRowMap(), AP->getColMap(), tpetra_prolong->getGlobalNumCols()));
  Teuchos::ArrayView<const global_ordinal_type> myLocalElements = BAP->getRowMap()->getLocalElementList();

  for (int color = 0; color < 3; ++color) {
    Teuchos::ArrayRCP<const scalar_type> localBAP = BAP_shrunk->getData(color);

    for (iter_type it = myLocalElements.begin(); it != myLocalElements.end(); ++it) {
      const local_ordinal_type i_local = *it;
      const local_ordinal_type aux     = BAP->getRowMap()->getLocalElement(i_local);

      std::vector<global_ordinal_type> BAP_inds;
      std::vector<scalar_type> BAP_vals;

      local_ordinal_type aux2;

      if ((mypid - 1) % 3 == color && (mypid - 1) >= 0 && (mypid - 1) < tpetra_prolong->getGlobalNumCols())
        aux2 = BAP->getColMap()->getLocalElement(mypid - 1);
      else if ((mypid - 2) % 3 == color && (mypid - 2) >= 0 && (mypid - 2) < tpetra_prolong->getGlobalNumCols())
        aux2 = BAP->getColMap()->getLocalElement(mypid - 2);
      else if ((mypid) % 3 == color && (mypid) >= 0 && (mypid) < tpetra_prolong->getGlobalNumCols())
        aux2 = BAP->getColMap()->getLocalElement(mypid);

      if (aux2 >= 0) {
        BAP_inds.emplace_back(aux2);
        BAP_vals.emplace_back(localBAP[aux]);
        BAP->insertLocalValues(aux, BAP_inds, BAP_vals);
      }
    }
  }
  BAP->fillComplete(tpetra_prolong->getDomainMap(), tpetra_prolong->getRangeMap());

  //=============================================================================================================

  RCP<crs_matrix_type> Pbar    = Tpetra::MatrixMatrix::add(1.0, false, *tpetra_prolong, -1.0, false, *BAP);
  RCP<xpetra_matrix> mueluPbar = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(Pbar);

  H->GetLevel(1)->Set("Pbar", mueluPbar);

  H->IsPreconditioner(true);
  B_coarse_ = rcp(new muelu_tpetra_operator_type(H));
}

void AdditiveVariant::apply(const multivector_type& r, multivector_type& Pr, Teuchos::ETransp mode = Teuchos::NO_TRANS, scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(), scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
  Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> Export_fine1(r.getMap(), DomainMap_);
  Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> Export_fine2(RangeMap_, r.getMap());

  multivector_type r_fine(DomainMap_, 1);
  r_fine.doImport(r, Export_fine1, Tpetra::INSERT);
  multivector_type B_fine_Pr(RangeMap_, 1);
  B_fine_->apply(r_fine, B_fine_Pr, mode, alpha, beta);
  multivector_type B_Pr1(r.getMap(), 1);
  B_Pr1.doImport(B_fine_Pr, Export_fine2, Tpetra::INSERT);

  multivector_type r_coarse(DomainMap_, 1);
  r_coarse.doImport(r, Export_fine1, Tpetra::INSERT);
  multivector_type B_coarse_Pr(RangeMap_, 1);
  B_coarse_->apply(r_coarse, B_coarse_Pr, mode, alpha, beta);
  multivector_type B_Pr2(r.getMap(), 1);
  B_Pr2.doImport(B_coarse_Pr, Export_fine2, Tpetra::INSERT);

  multivector_type B_Pr_sum(r.getMap(), 1);
  B_Pr_sum.update(1.0, B_Pr1, 1.0, B_Pr2, 0.0);  // Careful to set the correct coefficients!!!

  Pr = B_Pr_sum;
}

}  // namespace MueLu
