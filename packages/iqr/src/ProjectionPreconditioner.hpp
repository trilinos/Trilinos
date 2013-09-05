#ifndef PROJECTION_PRECONDITIONER_H
#define PROJECTION_PRECONDITIONER_H

#include <vector>

#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziConfigDefs.hpp>
#include <AnasaziEpetraAdapter.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include <Epetra_Import.h>

#include <EpetraExt_OperatorOut.h>
#include <EpetraExt_MultiVectorOut.h>

namespace IQR
{
/*
 * This class implements a preconditioner based on projection to a
 * subspace as follows:
 *
 * Let A x = b be the linear system that is solved with preconditioned
 * GMRES.
 *
 * G is an MxM matrix that is somewhat similar (?) to A (also MxM).
 *
 * We compute the N largest eigenvalues of G and their associated eigenvectors.
 *
 * We construct the matrix U1 (of size MxN) whose columns are the computed
 * eigenvectors.
 *
 * We construct a second MxN matrix V1 = A * U1
 *
 * The columns of V1 are orthonormalized with a Gram-Schmidt process resulting
 * in a new matrix V.
 * The projection coefficients and normalization of V1 are applied also to U1,
 * as follows:
 *
 *  V = V1;
 *  U = U1;
 *  for i = 1:s
 *      for j = 1:i-1
 *          p(i,j) = (V1(:,i)' * V(:,j)) / (V(:,j)' * V(:,j));
 *          V(:,i) = V(:,i) - (p(i,j) .* V(:,j));
 *          U(:,i) = U(:,i) - (p(i,j) .* U(:,j));
 *      end
 *  end
 *
 *  % Normalize
 *  for i = 1:s
 *      n = sqrt((V(:,i)' * V(:,i)));
 *      V(:,i) = V(:,i) ./ n;
 *      U(:,i) = U(:,i) ./ n;
 *  end
 *
 *  The columns of the matrix V represent the vectors of a basis of a subspace
 *  in which the linear system is solved
 *  The application of the preconditioner to a vector g is done in the following
 *  steps:
 *  1. The vector g is projected onto the subspace V and the inverse of the
 *  restriction to V of the matrix A is applied:
 *  A_V^{-1} * V * V' * g = U * V' * g
 *  2. The complementary projection of g is scaled by a value lambda which is
 *  chosen as the smallest eigenvalue computed, as an approximation of the mean
 *  eigenvalues of the matrix G
 *  1/lambda * (g - V * V' * g)
 *
 *  The full expression of the application of the preconditioner to a vector g
 *  reads:
 *
 *  gv = V' * g
 *
 *  gq = g - V * gv
 *
 *  P^{-1} g = U * gv + (1 / lambda) * gq;
 */
template < typename Map,
		   typename Operator,
		   typename MultiVector >
class ProjectionPreconditioner
{
public:
	// Public typedefs
	typedef MultiVector MV;
	typedef Operator OP;
	typedef Anasazi::MultiVecTraits<double, MV> MVT;
	typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

	// Constructor
	//
	ProjectionPreconditioner(Teuchos::RCP<const Operator> S,
							 Teuchos::RCP<const Operator> G,
							 const Teuchos::ParameterList parameters);
	// Destructor
	~ProjectionPreconditioner();

	// Public methods
	int Setup();
	int ApplyInverse(const MultiVector &B, MultiVector& X) const;

private:
	// Private methods
	int GramSchmidtCustom();

	// Data
	Teuchos::RCP<const Operator> S_;
	Teuchos::RCP<const Operator> G_;
	Teuchos::RCP<const Map> Smap_;
	Teuchos::RCP<const Map> Gmap_;
	Teuchos::RCP<MultiVector> Unonrestricted_;
	Teuchos::RCP<MultiVector> U_;
	Teuchos::RCP<MultiVector> V_;
	std::vector<Anasazi::Value<double> > eigenvalues_;
	Teuchos::ParameterList parameters_;
	double lambda_;
	int numVectors_;
	int myPID_;
};

template < typename Map,
		   typename Operator,
		   typename MultiVector >
ProjectionPreconditioner<Map, Operator, MultiVector>
		::ProjectionPreconditioner(Teuchos::RCP<const Operator> S,
			 	 	 	 	 	   Teuchos::RCP<const Operator> G,
			 	 	 	 	 	   const Teuchos::ParameterList parameters)
		 : S_(S),
		   G_(G),
		   Smap_(&(S_->OperatorDomainMap()), false),
		   Gmap_(&(G_->OperatorDomainMap()), false),
		   parameters_(parameters),
		   lambda_(1.0),
		   numVectors_(0),
		   myPID_(G_->Comm().MyPID())
{
}

template < typename Map,
		   typename Operator,
		   typename MultiVector >
ProjectionPreconditioner<Map, Operator, MultiVector>
		::~ProjectionPreconditioner()
{
}

template < typename Map,
		   typename Operator,
		   typename MultiVector >
int ProjectionPreconditioner<Map, Operator, MultiVector>
		::Setup()
{
	// Compute the largest N_ eigenvalues of G_ and the associated
	// eigenvectors: [U1, D] = eigs(G_, N_)

	double relativeSubspaceDim = parameters_.get<double>("relative subspace dimension", 0.5);
	numVectors_ = static_cast<int>(relativeSubspaceDim * Smap_->NumGlobalElements());

	Teuchos::RCP<MultiVector> ivec = Teuchos::rcp(new MultiVector(*Gmap_, 1));
	ivec->Random();

	Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem;
	MyProblem = Teuchos::rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(G_, ivec));
	MyProblem->setHermitian(false);
	MyProblem->setNEV(numVectors_);

	bool status = MyProblem->setProblem();
	if (status != true) {
		if (myPID_ == 0) {
			cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
		}
	#ifdef HAVE_MPI
		MPI_Finalize() ;
	#endif
		return -1;
	}

	Teuchos::ParameterList& AnasaziPL = parameters_.sublist("anasazi parameters", false);
	Teuchos::ParameterList tempList;
	tempList.set( "Verbosity", Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary);
	tempList.set( "Which", "LM" );
	tempList.set( "Maximum Restarts", 500 );
	tempList.set( "Step Size", 5 );
	AnasaziPL.setParametersNotAlreadySet(tempList);

	Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, AnasaziPL);

	Anasazi::ReturnType returnCode = MySolverMgr.solve();
	if (returnCode != Anasazi::Converged && myPID_ == 0) {
		cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
	}

	Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
	eigenvalues_ = sol.Evals;
	Unonrestricted_ = sol.Evecs;
	std::vector<int> index = sol.index;

	if (sol.numVecs != numVectors_) {
		if (!myPID_) {
			std::cout << "Number of computed eigenvalues lower than requested" << std::endl;
		}
#ifdef HAVE_MPI
		MPI_Finalize() ;
#endif
		return -1;
	}

	// If S_ and G_ don't have the same map we need to restrict the U_ vectors
	// to the same space as S_
	if (! Smap_->SameAs(*Gmap_)) {
		// We need to build an importer, restric U_ to the map of S_ and
		Epetra_Import restrictor(*Smap_, *Gmap_);
		U_ = Teuchos::rcp(new MultiVector(*Smap_, numVectors_));
		U_->Import(*Unonrestricted_, restrictor, Insert);
	} else {
		// U_ and Unonrestricted_ should point to the same thing
		U_ = Unonrestricted_;
	}

    // Compute V1 = A_ * U1
	V_ = Teuchos::rcp(new MultiVector(*Smap_, numVectors_));
    for (int i = 0; i < numVectors_; ++i) {
        MultiVector& ui = *((*U_)(i));
        MultiVector& vi = *((*V_)(i));
        S_->Apply(ui, vi);
    }

	// Compute [V, U] = GramSchmidtCustom(V1, U1)
	GramSchmidtCustom();

	// Compute lambda = D(N_ - 1)
	lambda_ = eigenvalues_[numVectors_ - 1].realpart;

	return 0;
}

template < typename Map,
		   typename Operator,
		   typename MultiVector >
int ProjectionPreconditioner<Map, Operator, MultiVector>
		::ApplyInverse(const MultiVector &B, MultiVector& X) const
{
    // TODO: you can actually get rid of the second loop in i
    // and also of the Bq temp vector
	for (int k = 0; k < B.NumVectors(); ++k) {
        const MultiVector* currentB = B(k);
		MultiVector* currentX = X(k);

		// Compute Bv = V^T * B
        // Compute Bq = (1 / lambda) * (B - V * Bv)
        std::vector<double> Bv(numVectors_, 0.0);
        Teuchos::RCP<MultiVector> Bq;
        Bq = MVT::CloneCopy(*currentB);
        for (int i = 0; i < numVectors_; ++i) {
            MultiVector& vi = *((*V_)(i));
            vi.Dot(*currentB, &Bv[i]);
            Bq->Update(-1.0 * Bv[i], vi, 1.0);
        }
        Bq->Scale(1.0 / lambda_);

		// X = U * Bv + Bq
		*currentX = *Bq;
		for (int i = 0; i < numVectors_; ++i) {
			MultiVector& ui = *((*U_)(i));
			currentX->Update(Bv[i], ui, 1.0);
		}
	}
	return 0;
}

template < typename Map,
		   typename Operator,
		   typename MultiVector >
int ProjectionPreconditioner<Map, Operator, MultiVector>::GramSchmidtCustom()
{
	// First the orthogonalization part
	std::vector<double> p(numVectors_ * numVectors_, 0.0);
	for (int i = 0; i < numVectors_; ++i) {
        MultiVector* vi = (*V_)(i);
        MultiVector* ui = (*U_)(i);
		for (int j = 0; j < i; ++j) {
			double& pij = p[i * numVectors_ + j];
			double a, b;
            MultiVector* vj = (*V_)(j);
            MultiVector* uj = (*U_)(j);
            vi->Dot(*vj, &a);
			vj->Dot(*vj, &b);
			pij = a / b;
			vi->Update(-1.0 * pij, *vj, 1.0);
			ui->Update(-1.0 * pij, *uj, 1.0);
		}
	}

	// Then the normalization
	std::vector<double> normsV(numVectors_, 0.0);
    MVT::MvNorm(*V_, normsV);
	for (int i = 0; i < numVectors_; ++i) {
        (*V_)(i)->Scale(1.0 / normsV[i]);
        (*U_)(i)->Scale(1.0 / normsV[i]);
	}

	return 0;
}

}

#endif // PROJECTION_PRECONDITIONER
