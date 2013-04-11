// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef GALERI_HELMHOLTZFEM2DPROBLEM_HPP
#define GALERI_HELMHOLTZFEM2DPROBLEM_HPP

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_Problem.hpp"
#include "Galeri_MultiVectorTraits.hpp"
#include "Galeri_XpetraUtils.hpp"

// Finite element discretization for the indefinite Helmholtz equation

namespace Galeri {

  namespace Xpetra {

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class HelmholtzFEM2DProblem : public Problem<Map,Matrix,MultiVector> {

    public:

      HelmholtzFEM2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) {

	// parameters
	hx_          = list.get("stretchx", 1.0);
	hy_          = list.get("stretchy", 1.0);
	px_          = list.get("px", 1);
	py_          = list.get("py", 1);
        nx_          = list.get("nx", -1);
        ny_          = list.get("ny", -1);
	omega_       = list.get("omega", 2.0*M_PI);
	shift_       = list.get("shift", 0.0);
	delta_       = list.get("delta", 2.0);
	PMLx_left    = list.get("PMLx_left",  0);
	PMLx_right   = list.get("PMLx_right", 0);
	PMLy_left    = list.get("PMLy_left",  0);
	PMLy_right   = list.get("PMLy_right", 0);	
	
	// calculate info
	Dx_    = hx_*nx_;
	Dy_    = hy_*ny_;
	LBx_   = PMLx_left*hx_;
	RBx_   = Dx_-PMLx_right*hx_;
	LBy_   = PMLy_left*hy_;
	RBy_   = Dy_-PMLy_right*hy_;
	PMLwidthx_ = max(LBx_,Dx_-RBx_); if(PMLwidthx_==0) {  PMLwidthx_=1.0; }
	PMLwidthy_ = max(LBy_,Dy_-RBy_); if(PMLwidthy_==0) {  PMLwidthy_=1.0; }
        nDim = 2;

	// currently AMG doesn't suport higher orders. set to linears
	px_=1;
	py_=1;

        // NOTE: galeri counts points, not elements
        dims.push_back(nx_*px_-1);
        dims.push_back(ny_*py_-1);
        TEUCHOS_TEST_FOR_EXCEPTION(nx_ <= 0 || ny_ <= 0, std::logic_error, "nx and ny must be positive");

      }

      Teuchos::RCP<Matrix>      BuildMatrix();

    private:

      typedef Scalar        SC;
      typedef LocalOrdinal  LO;
      typedef GlobalOrdinal GO;

      struct Point {
        double x, y, z;
        Point() { z = 0.0; }
        Point(double x_, double y_, double z_ = 0.0) : x(x_), y(y_), z(z_) { }
      };

      // General mesh parameters
      GlobalOrdinal                   nx_, ny_, nz_;
      size_t                          nDim, px_, py_, pz_;
      std::vector<GO>                 dims;
      std::vector<Point>              nodes;
      std::vector< std::vector<LO> >  elements;
      std::vector<GO>                 local2Global_;

      // Helmholtz/PML parameters
      double   hx_, hy_, hz_, shift_, delta_, omega_;
      int      PMLx_left, PMLx_right; 
      int      PMLy_left, PMLy_right;
      int      PMLz_left, PMLz_right;
      double   Dx_, Dy_, Dz_;
      double   LBx_, RBx_, LBy_, RBy_, LBz_, RBz_;
      double   PMLwidthx_, PMLwidthy_, PMLwidthz_;

      void BuildMesh();
      void BuildPoints(std::vector<Point>& quadPoints, std::vector<double>& quadWeights);
      void EvalBasis(Point& quadPoint, std::vector<double>& vals, std::vector<double>& dxs, std::vector<double>& dys);
      void EvalStretch(double& shiftx, double& shifty, const std::vector<Point>& quadPoints, std::vector<Scalar>& sx, std::vector<Scalar>& sy);

    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> HelmholtzFEM2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      using Teuchos::SerialDenseMatrix;

      BuildMesh();

      const size_t numDofPerNode   = 1;
      const size_t numNodesPerElem = (px_+1)*(py_+1);
      const size_t numDofPerElem   = numNodesPerElem * numDofPerNode;

      TEUCHOS_TEST_FOR_EXCEPTION(elements[0].size() != numNodesPerElem, std::logic_error, "Incorrect number of element vertices");

      // Compute quadrature points/rules
      std::vector<Point>  quadPoints;
      std::vector<double> quadWeights;
      BuildPoints(quadPoints, quadWeights);

      // Compute basis function values and derivatives at each quadrature point
      std::vector< vector<double> > vals, dxs, dys;
      vals.resize(quadPoints.size());
      dxs.resize(quadPoints.size());
      dys.resize(quadPoints.size());
      for(int i=0; i<quadPoints.size(); i++) {
	EvalBasis(quadPoints[i], vals[i], dxs[i], dys[i]);
      }

      this->A_ = MatrixTraits<Map,Matrix>::Build(this->Map_, (2*px_+1)*(2*py_+1)*numDofPerNode);
      SC one = Teuchos::ScalarTraits<SC>::one();
      SC zero = Teuchos::ScalarTraits<SC>::zero();
      SC cpxshift(1.0,shift_);

      // iterate over elements
      for (size_t i = 0; i < elements.size(); i++) {
	
        SerialDenseMatrix<LO,SC> KE(numDofPerElem, numDofPerElem);

	// element domain is [shiftx,shiftx+hx] x [shifty,shifty+hy]
        std::vector<LO>& elemNodes = elements[i];
	double shiftx = nodes[elemNodes[0]].x;
	double shifty = nodes[elemNodes[0]].y;

	// evaluate PML values at quadrature points for this element
	std::vector<SC> stretchx, stretchy;
	EvalStretch(shiftx, shifty, quadPoints, stretchx, stretchy);

        // Evaluate the stiffness matrix for the element
        for (size_t j = 0; j < quadPoints.size(); j++) {
	  Scalar sx=stretchx[j];
	  Scalar sy=stretchy[j];
	  double qdwt=quadWeights[j];
	  std::vector<double> curvals=vals[j];
	  std::vector<double> curdxs=dxs[j];
	  std::vector<double> curdys=dys[j];
	  Scalar mass = qdwt*cpxshift*omega_*omega_*sx*sy;
	  Scalar pml1 = qdwt*sy/sx;
	  Scalar pml2 = qdwt*sx/sy;
	  for(int m=0; m<numDofPerElem; m++) {
	    for(int n=0; n<numDofPerElem; n++) {
	      KE[m][n] += pml1*curdxs[m]*curdxs[n] + pml2*curdys[m]*curdys[n] - mass*curvals[m]*curvals[n];
	    }
	  }
        }

	Teuchos::Array<GO> elemDofs(numDofPerElem);
	for (size_t j = 0; j < numDofPerElem; j++) {
	  elemDofs[j] = local2Global_[elemNodes[j]];
	  if(local2Global_[elemNodes[j]]>9999) {
	    std::cout<<"index "<<local2Global_[elemNodes[j]]<<std::endl;
	  }	    
	}

        // Insert KE into the global matrix
        for (size_t j = 0; j < numDofPerElem; j++)
          if (this->Map_->isNodeGlobalElement(elemDofs[j]))
            this->A_->insertGlobalValues(elemDofs[j], elemDofs, Teuchos::ArrayView<SC>(KE[j], numDofPerElem));

      }

      this->A_->fillComplete();
      return this->A_;

    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMesh() {

      GO myPID  = this->Map_->getComm()->getRank();
      GO mySize = this->Map_->getComm()->getSize();
      GO nx = -1,                       ny = -1;
      GO mx = this->list_.get("mx", 1), my = this->list_.get("my", 1);
      if(mx*my != mySize) {
	int squareRoot = std::max(1,(int)std::floor(sqrt((double) mySize)));
	my = squareRoot;
	while( mySize % my != 0 )
	  ++my;
	mx = mySize / my;
      }
      GO shiftx, shifty;

      Utils::getSubdomainData(dims[0], mx, myPID % mx, nx, shiftx);
      Utils::getSubdomainData(dims[1], my, myPID / mx, ny, shifty);

      // Expand subdomain to do overlap
      if (shiftx    > 0)        { nx++; shiftx--; }
      if (shifty    > 0)        { ny++; shifty--; }
      if (shiftx+nx < dims[0])  { nx++;           }
      if (shifty+ny < dims[1])  { ny++;           }

      nodes        .resize((nx+1)*(ny+1));
      local2Global_.resize((nx+1)*(ny+1));
      elements     .resize(nx*ny);

#define NODE(i,j) ((j)*(nx+1) + (i))
#define CELL(i,j) ((j)*nx     + (i))
      for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) {
          int ii = shiftx + i, jj = shifty + j;
          int nodeID = NODE(i,j);
          nodes[nodeID] = Point((ii+1)*hx_, (jj+1)*hy_);
          local2Global_[nodeID] = jj*nx_ + ii;
        }

      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
          std::vector<int>& element = elements[CELL(i,j)];
          element.resize(4);
          element[0] = NODE(i,  j  );
          element[1] = NODE(i+1,j  );
          element[2] = NODE(i,  j+1);
          element[3] = NODE(i+1,j+1);
        }
#undef NODE
#undef CELL

    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildPoints(std::vector<Point>& quadPoints, std::vector<double>& quadWeights) {

      // Use nonstandard quadrature for Helmholtz from the following paper:
      // M. Ainsworth, H.A. Wajid. Optimally blended spectral-finite element scheme for wave propagation
      // and nonstandard reduced integration. SIAM J. Numer. Analysis, 48(1): 346-371, 2010.

      // Make reference domain [0,hx_] x [0,hy_].
      double qdpt = 0.816496580927726;
      quadPoints.resize(4);
      quadPoints[0] = Point( (1.0+qdpt)*hx_/2.0, (1.0+qdpt)*hy_/2.0 );
      quadPoints[1] = Point( (1.0+qdpt)*hx_/2.0, (1.0-qdpt)*hy_/2.0 );
      quadPoints[2] = Point( (1.0-qdpt)*hx_/2.0, (1.0+qdpt)*hy_/2.0 );
      quadPoints[3] = Point( (1.0-qdpt)*hx_/2.0, (1.0-qdpt)*hy_/2.0 );
      quadWeights.resize(4);
      quadWeights[0] = hx_*hy_/4.0;
      quadWeights[1] = hx_*hy_/4.0;
      quadWeights[2] = hx_*hy_/4.0;
      quadWeights[3] = hx_*hy_/4.0;

    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalBasis(Point& quadPoint, std::vector<double>& vals, std::vector<double>& dxs, std::vector<double>& dys) {
      // For the reference element, evaluate each basis function and its derivatives at a particular quadrature point
      int numBasis=(px_+1)*(py_+1);
      vals.resize(numBasis);
      dxs.resize(numBasis);
      dys.resize(numBasis);
      double x=quadPoint.x;
      double y=quadPoint.y;
      double valx1=(hx_-x)/hx_;   double dx1=-1/hx_;
      double valx2=x/hx_;         double dx2=1/hx_;
      double valy1=(hy_-y)/hy_;   double dy1=-1/hy_;
      double valy2=y/hy_;         double dy2=1/hy_;
      // linear function values, derivatives in x, y directions
      vals[0]=valx1*valy1;  dxs[0]=dx1*valy1;  dys[0]=valx1*dy1;
      vals[1]=valx2*valy1;  dxs[1]=dx2*valy1;  dys[1]=valx2*dy1;
      vals[2]=valx1*valy2;  dxs[2]=dx1*valy2;  dys[2]=valx1*dy2;
      vals[3]=valx2*valy2;  dxs[3]=dx2*valy2;  dys[3]=valx2*dy2;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalStretch(double& shiftx, double& shifty, const std::vector<Point>& quadPoints, std::vector<Scalar>& sx, std::vector<Scalar>& sy) {
      // For the current element domain, evaluate the PML stretching functions at the quadrature points
      sx.resize(quadPoints.size());
      sy.resize(quadPoints.size());
      for(int i=0; i<quadPoints.size(); i++) {
	double quadx=quadPoints[i].x;
	double quady=quadPoints[i].y;
	double curx=shiftx+quadx;
	double cury=shifty+quady;
	double sigx, sigy;
	if(curx<LBx_)        { sigx = delta_*pow((curx-LBx_)/PMLwidthx_,2);  }
	else if(curx>RBx_)   { sigx = delta_*pow((curx-RBx_)/PMLwidthx_,2);  }
	else                 { sigx = 0.0;                                   }
	if(cury<LBy_)        { sigy = delta_*pow((cury-LBy_)/PMLwidthy_,2);  }
	else if(cury>RBy_)   { sigy = delta_*pow((cury-RBy_)/PMLwidthy_,2);  }
	else                 { sigy = 0.0;                                   }
	Scalar pmlx(1.0,sigx);
	Scalar pmly(1.0,sigy);
	sx[i]=pmlx;
	sy[i]=pmly;
      }
    }


  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_HELMHOLTZFEM2DPROBLEM_HPP
