// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_PDEOPT_ELASTICITYSIMP_OPERATORS_H
#define ROL_PDEOPT_ELASTICITYSIMP_OPERATORS_H

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/elasticitySIMP.hpp"

template<class Real>
class ElasticitySIMPOperators : public ElasticitySIMP <Real> {

private:

Real volFrac_;

public:

ElasticitySIMPOperators(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
             const Teuchos::RCP<Teuchos::ParameterList> &parlist,
             const Teuchos::RCP<std::ostream> &outStream) 
{
	this->ElasticitySIMPOperators_Initialize (comm, parlist, outStream);
	this->SetSIMPParallelStructure();
	this->SetUpLocalIntrepidArrays();
	this->ComputeLocalSystemMats(true);
	this->ComputeLocalForceVec();
	std::vector<int> dbc_side {0};
	this->SetUpMyDBCInfo(dbc_side);
	
	this->AssembleSystemMats();
	this->AssembleRHSVector();
	
	this->EnforceDBC();
	this->ConstructSolvers();
}

void ElasticitySIMPOperators_Initialize(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
        	   	const Teuchos::RCP<Teuchos::ParameterList> &parlist,
        	   	const Teuchos::RCP<std::ostream> &outStream) 
{
	this->ElasticitySIMP_Initialize(comm, parlist, outStream);
	volFrac_     = this->parlist_->sublist("ElasticityTopoOpt").get("volfrac", 0.5);
	this->initDensity_ = volFrac_;
}

// construct solvers with new material
void constructSolverWithNewMaterial()
{
	bool ifInit = false;
	this->ComputeLocalSystemMats(ifInit);
	this->AssembleSystemMats();
	this->EnforceDBC();
	this->ConstructSolvers();
}
//

Teuchos::RCP<Intrepid::FieldContainer<int> > getPosCell() const 
{
    	return this->meshMgr_->getPosCell();
}

void ApplyMatAToVec (Teuchos::RCP<Tpetra::MultiVector<> > &Ju, const Teuchos::RCP<const Tpetra::MultiVector<> > &u)
{
	this->matA_dirichlet_->apply(*u, *Ju);
}

void ApplyJacobian1ToVec (Teuchos::RCP<Tpetra::MultiVector<> > &Ju, const Teuchos::RCP<const Tpetra::MultiVector<> > &u)
{
	// Jv is myUniqueMap_
	// u is myUniqueMap_
	// v should be myCellMap_
	//
	Teuchos::RCP<Tpetra::MultiVector<> > Ju_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Teuchos::RCP<Tpetra::MultiVector<> > u_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Tpetra::Import<> importer(u->getMap(), u_overlap->getMap());
	u_overlap->doImport(*u, importer, Tpetra::REPLACE);
        //simply apply matrix to a vector, no need to apply EBC to the vec
	//ApplyBCToVec (u_overlap);
	Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();
	
	Real SIMPScale;
	Real sum;
	Intrepid::FieldContainer<Real> localDisp(this->numLocalDofs_);

	for(int i=0; i<this->numCells_; i++)
	{
	    SIMPScale = this->SIMPmaterial_[i]->getSIMPScaleFactor();
	    localDisp.initialize(0.0);
	    for(int j=0; j<this->numLocalDofs_; j++)
		localDisp(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];

	    for(int j=0; j<this->numLocalDofs_; j++)
	    {	
		if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j)))
			Ju_overlap -> replaceGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, localDisp(j));	
		else
		{
			sum = 0.0;
			for(int k=0; k<this->numLocalDofs_; k++)
			{
				if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) )
				{
    					sum += SIMPScale * (*this->gradgradMats0_)(i, j, k) * localDisp(k);
				}
			}
			Ju_overlap -> sumIntoGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, sum);	
		}
	    }
	}
	
	Tpetra::Export<> exporter(Ju_overlap->getMap(), Ju->getMap());
	Ju->doExport(*Ju_overlap, exporter, Tpetra::REPLACE);
}

void ApplyInverseJacobian1ToVec (Teuchos::RCP<Tpetra::MultiVector<> > &InvJu, const Teuchos::RCP<const Tpetra::MultiVector<> > &u, bool ifTrans)
{
  	//Ju is matA->getDomainMap();
    	this->getSolver(ifTrans)->setB(u);
    	this->getSolver(ifTrans)->setX(InvJu);
    	this->getSolver(ifTrans)->solve();
}

void ApplyJacobian2ToVec (Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
			const Teuchos::RCP<const Tpetra::MultiVector<> > &u, const Teuchos::RCP<const Tpetra::MultiVector<> > &v)
{
	// Jv is myUniqueMap_
	// u is myUniqueMap_
	// v should be myCellMap_
	//
	Teuchos::RCP<Tpetra::MultiVector<> > Jv_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Teuchos::RCP<Tpetra::MultiVector<> > u_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Tpetra::Import<> importer(u->getMap(), u_overlap->getMap());
	u_overlap->doImport(*u, importer, Tpetra::REPLACE);
	//apply BC here, KU
        this->ApplyBCToVec (u_overlap);
	Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();
	
	Teuchos::RCP<Tpetra::MultiVector<> > v_local = Tpetra::rcp(new Tpetra::MultiVector<>(this->myCellMap_, 1, true));
	Tpetra::Export<> exporter1(v->getMap(), this->myCellMap_);
	v_local->doExport(*v, exporter1, Tpetra::REPLACE);
	Teuchos::ArrayRCP<const Real> vData = v_local->get1dView();
	
	Real SIMPDScale;
	Real localV;
	Real sum;
	Intrepid::FieldContainer<Real> localDisp(this->numLocalDofs_);

	for(int i=0; i<this->numCells_; i++)
	{
	    SIMPDScale = this->SIMPmaterial_[i]->getSIMPFirstDerivativeScaleFactor();
	    localV = vData[v_local->getMap()->getLocalElement(this->myCellIds_[i])];
	
	    localDisp.initialize(0.0);
	    for(int j=0; j<this->numLocalDofs_; j++)
		localDisp(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];

	    for(int j=0; j<this->numLocalDofs_; j++)
	    {	
		if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j)))
			Jv_overlap -> replaceGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, localDisp(j));	
		else
		{
			sum = 0.0;
			for(int k=0; k<this->numLocalDofs_; k++)
			{
				if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) )
    					sum += localV * SIMPDScale * (*this->gradgradMats0_)(i, j, k) * localDisp(k);
			}
			
			Jv_overlap -> sumIntoGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, sum);
		}
	    }
	}
	
	Tpetra::Export<> exporter2(Jv_overlap->getMap(), Jv->getMap());
	Jv->doExport(*Jv_overlap, exporter2, Tpetra::REPLACE);
}

	
void ApplyAdjointJacobian2ToVec (Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
			const Teuchos::RCP<const Tpetra::MultiVector<> > &u, const Teuchos::RCP<const Tpetra::MultiVector<> > &v)
{
	// Jv is myCellMap_
	// u should be myUniqueMap_
	// v should be myUniqueMap_
	//
	Teuchos::RCP<Tpetra::MultiVector<> > u_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Tpetra::Import<> importer1(u->getMap(), u_overlap->getMap());
	u_overlap->doImport(*u, importer1, Tpetra::REPLACE);
	//only apply BC to u
	this->ApplyBCToVec (u_overlap);
	Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();
	
	Teuchos::RCP<Tpetra::MultiVector<> > v_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Tpetra::Import<> importer2(v->getMap(), v_overlap->getMap());
	v_overlap->doImport(*v, importer2, Tpetra::REPLACE);
	//ApplyBCToVec (v_overlap);
	Teuchos::ArrayRCP<const Real> vData = v_overlap->get1dView();
    	
	Intrepid::FieldContainer<Real> u_local(this->numLocalDofs_);
	Intrepid::FieldContainer<Real> v_local(this->numLocalDofs_);
	
	Real SIMPDScale;
	Real sum;
	
	for(int i=0; i<this->numCells_; i++)
	{
	    SIMPDScale = this->SIMPmaterial_[i]->getSIMPFirstDerivativeScaleFactor();
	    u_local.initialize(0.0);
	    v_local.initialize(0.0);
	    
	    for(int j=0; j<this->numLocalDofs_; j++)
	    {
		u_local(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
		v_local(j) = vData[v_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
	    }

	    sum = 0.0;	
	    for(int j=0; j<this->numLocalDofs_; j++)
	    {	
		if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j)))
			sum += u_local(j) * v_local(j);
		else
		{
			for(int k=0; k<this->numLocalDofs_; k++)
			{
				if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) )
    					sum += SIMPDScale * (*this->gradgradMats0_)(i, j, k) * v_local(j) * u_local(k);
			}
		}
	    }
	    //put value into myDensityDerivative_ 
	    Jv -> replaceGlobalValue(this->myCellIds_[i], 0, sum);	
	}
}

void ApplyAdjointHessian22ToVec (Teuchos::RCP<Tpetra::MultiVector<> > &Hvw, const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
			const Teuchos::RCP<const Tpetra::MultiVector<> > &v, const Teuchos::RCP<const Tpetra::MultiVector<> > &w)
{
	// Hvw is myCellMap_
	// u should be myUniqueMap_
	// v should be myCellMap_
	// w shluld be myUniqueMapi_
	//
	Teuchos::RCP<Tpetra::MultiVector<> > u_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Tpetra::Import<> importer1(u->getMap(), u_overlap->getMap());
	u_overlap->doImport(*u, importer1, Tpetra::REPLACE);
	//only apply BC to U
	this->ApplyBCToVec (u_overlap);
	Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();
	
	Tpetra::Export<> exporter(v->getMap(), this->myCellMap_);
	Teuchos::RCP<Tpetra::MultiVector<> > v_local = Tpetra::rcp(new Tpetra::MultiVector<>(this->myCellMap_, 1, true));
	v_local->doExport(*v, exporter, Tpetra::REPLACE);
	Teuchos::ArrayRCP<const Real> vData = v_local->get1dView();
	
	Teuchos::RCP<Tpetra::MultiVector<> > w_overlap = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
	Tpetra::Import<> importer2(w->getMap(), w_overlap->getMap());
	w_overlap->doImport(*w, importer2, Tpetra::REPLACE);
	//ApplyBCToVec (w_overlap);
	Teuchos::ArrayRCP<const Real> wData = w_overlap->get1dView();

	Intrepid::FieldContainer<Real> u_local(this->numLocalDofs_);
	Intrepid::FieldContainer<Real> w_local(this->numLocalDofs_);
	
	Real SIMPDDScale;
	Real localV;
	Real sum;
	
	for(int i=0; i<this->numCells_; i++)
	{
	    SIMPDDScale = this->SIMPmaterial_[i]->getSIMPSecondDerivativeScaleFactor();
	    localV = vData[v_local->getMap()->getLocalElement(this->myCellIds_[i])];
	
	    u_local.initialize(0.0);
	    w_local.initialize(0.0);
	    
	    for(int j=0; j<this->numLocalDofs_; j++)
	    {
		u_local(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
		w_local(j) = wData[w_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
	    }

	    sum = 0.0;	
	    for(int j=0; j<this->numLocalDofs_; j++)
	    {	
		if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j)))
			sum += u_local(j) * w_local(j);
		else
		{
			for(int k=0; k<this->numLocalDofs_; k++)
			{
				if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) )
    					sum += localV * SIMPDDScale * (*this->gradgradMats0_)(i, j, k) * w_local(j) * u_local(k);
			}
		}
	    }
	    //put value into myDensityDerivative_ 
	    Hvw -> replaceGlobalValue(this->myCellIds_[i], 0, sum);	
	}
}

}; // class ElasticitySIMPOperators

#endif
