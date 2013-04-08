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

#ifndef GALERI_HELMHOLTZFEM3DPROBLEM_HPP
#define GALERI_HELMHOLTZFEM3DPROBLEM_HPP

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_Problem.hpp"
#include "Galeri_MultiVectorTraits.hpp"
#include "Galeri_XpetraUtils.hpp"

// Finite element discretization for the indefinite Helmholtz equation

namespace Galeri {

  namespace Xpetra {

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class HelmholtzFEM3DProblem : public Problem<Map,Matrix,MultiVector> {

    public:

      HelmholtzFEM3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) {

	// parameters
	hx_          = list.get("stretchx", 1.0);
	hy_          = list.get("stretchy", 1.0);
	hz_          = list.get("stretchz", 1.0);
	px_          = list.get("px", 1);
	py_          = list.get("py", 1);
	pz_          = list.get("pz", 1);
        nx_          = list.get("nx", -1);
        ny_          = list.get("ny", -1);
        nz_          = list.get("nz", -1);
	omega_       = list.get("omega", 2.0*M_PI);
	shift_       = list.get("shift", 0.0);
	delta_       = list.get("delta", 2.0);
	PMLx_left    = list.get("PMLx_left",  0);
	PMLx_right   = list.get("PMLx_right", 0);
	PMLy_left    = list.get("PMLy_left",  0);
	PMLy_right   = list.get("PMLy_right", 0);	
	PMLz_left    = list.get("PMLz_left",  0);
	PMLz_right   = list.get("PMLz_right", 0);	
	
	// calculate info
	Dx_    = hx_*nx_;
	Dy_    = hy_*ny_;
	Dz_    = hz_*nz_;
	LBx_   = PMLx_left*hx_;
	RBx_   = Dx_-PMLx_right*hx_;
	LBy_   = PMLy_left*hy_;
	RBy_   = Dy_-PMLy_right*hy_;
	LBz_   = PMLz_left*hz_;
	RBz_   = Dz_-PMLz_right*hz_;
	PMLwidthx_ = max(LBx_,Dx_-RBx_); if(PMLwidthx_==0) {  PMLwidthx_=1.0; }
	PMLwidthy_ = max(LBy_,Dy_-RBy_); if(PMLwidthy_==0) {  PMLwidthy_=1.0; }
	PMLwidthz_ = max(LBz_,Dz_-RBz_); if(PMLwidthz_==0) {  PMLwidthz_=1.0; }
        nDim = 3;

	// currently AMG doesn't suport higher orders. set to linears
	px_=1;
	py_=1;
	pz_=1;

        // NOTE: galeri counts points, not elements
        dims.push_back(nx_*px_-1);
        dims.push_back(ny_*py_-1);
        dims.push_back(nz_*pz_-1);
        TEUCHOS_TEST_FOR_EXCEPTION(nx_ <= 0 || ny_ <= 0 || nz_ <= 0, std::logic_error, "nx, ny, and nz must be positive");

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
      void EvalBasis(Point& quadPoint, std::vector<double>& vals, std::vector<double>& dxs, std::vector<double>& dys, std::vector<double>& dzs);
      void EvalStretch(double& shiftx, double& shifty, double& shiftz, const std::vector<Point>& quadPoints, std::vector<Scalar>& sx, std::vector<Scalar>& sy, std::vector<Scalar>& sz);

    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> HelmholtzFEM3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      using Teuchos::SerialDenseMatrix;

      BuildMesh();

      const size_t numDofPerNode   = 1;
      const size_t numNodesPerElem = (px_+1)*(py_+1)*(pz_+1);
      const size_t numDofPerElem   = numNodesPerElem * numDofPerNode;

      TEUCHOS_TEST_FOR_EXCEPTION(elements[0].size() != numNodesPerElem, std::logic_error, "Incorrect number of element vertices");

      // Compute quadrature points/rules
      std::vector<Point>  quadPoints;
      std::vector<double> quadWeights;
      BuildPoints(quadPoints, quadWeights);

      // Compute basis function values and derivatives at each quadrature point
      std::vector< vector<double> > vals, dxs, dys, dzs;
      vals.resize(quadPoints.size());
      dxs.resize(quadPoints.size());
      dys.resize(quadPoints.size());
      dzs.resize(quadPoints.size());
      for(int i=0; i<quadPoints.size(); i++) {
	EvalBasis(quadPoints[i], vals[i], dxs[i], dys[i], dzs[i]);
      }

      this->A_ = MatrixTraits<Map,Matrix>::Build(this->Map_, (2*px_+1)*(2*py_+1)*(2*pz_+1)*numDofPerNode);
      SC one = Teuchos::ScalarTraits<SC>::one();
      SC zero = Teuchos::ScalarTraits<SC>::zero();
      SC cpxshift(1.0,shift_);

      // iterate over elements
      for (size_t i = 0; i < elements.size(); i++) {
	
        SerialDenseMatrix<LO,SC> KE(numDofPerElem, numDofPerElem);

	// element domain is [shiftx,shiftx+hx] x [shifty,shifty+hy] x [shiftz,shiftz+hz]
        std::vector<LO>& elemNodes = elements[i];
	double shiftx = nodes[elemNodes[0]].x;
	double shifty = nodes[elemNodes[0]].y;
	double shiftz = nodes[elemNodes[0]].z;

	// evaluate PML values at quadrature points for this element
	std::vector<SC> stretchx, stretchy, stretchz;
	EvalStretch(shiftx, shifty, shiftz, quadPoints, stretchx, stretchy, stretchz);

        // Evaluate the stiffness matrix for the element
        for (size_t j = 0; j < quadPoints.size(); j++) {
	  Scalar sx=stretchx[j];
	  Scalar sy=stretchy[j];
	  Scalar sz=stretchz[j];
	  double qdwt=quadWeights[j];
	  std::vector<double> curvals=vals[j];
	  std::vector<double> curdxs=dxs[j];
	  std::vector<double> curdys=dys[j];
	  std::vector<double> curdzs=dzs[j];
	  Scalar mass = qdwt*cpxshift*omega_*omega_*sx*sy*sz;
	  Scalar pml1 = qdwt*sy*sz/sx;
	  Scalar pml2 = qdwt*sx*sz/sy;
	  Scalar pml3 = qdwt*sx*sy/sz;
	  for(int m=0; m<numDofPerElem; m++) {
	    for(int n=0; n<numDofPerElem; n++) {
	      KE[m][n] += pml1*curdxs[m]*curdxs[n] + pml2*curdys[m]*curdys[n] + pml3*curdzs[m]*curdzs[n] - mass*curvals[m]*curvals[n];
	    }
	  }
        }

	Teuchos::Array<GO> elemDofs(numDofPerElem);
	for (size_t j = 0; j < numDofPerElem; j++) {
	  elemDofs[j] = local2Global_[elemNodes[j]];
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
    void HelmholtzFEM3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMesh() {

      GO myPID  = this->Map_->getComm()->getRank();
      GO mySize = this->Map_->getComm()->getSize();
      GO nx = -1,                        ny = -1,                        nz = -1;
      GO mx = this->list_.get("mx", 1),  my = this->list_.get("my", 1),  mz = this->list_.get("mz", 1);
      if(mx*my*mz != mySize) {
	const int cubeRoot = std::max(1,(int)std::floor(pow((double) mySize, 1./3.)));
	mz = cubeRoot;
	while( mySize % mz != 0 )
	  ++mz;
	const int xyCommSize = mySize / mz;
	const int squareRoot = std::max(1,(int)std::floor(sqrt((double) xyCommSize)));
	my = squareRoot;
	while( xyCommSize % my != 0 )
	  ++my;
	mx = xyCommSize / my;
      }
      GO shiftx, shifty, shiftz;

      Utils::getSubdomainData(dims[0], mx, myPID % mx, nx, shiftx);
      Utils::getSubdomainData(dims[1], my, ((myPID - (mx*my) * (myPID / (mx*my))) / mx), ny, shifty);
      Utils::getSubdomainData(dims[2], mz, myPID / (mx*my), nz, shiftz);

      // Expand subdomain to do overlap
      if (shiftx    > 0)        { nx++; shiftx--; }
      if (shifty    > 0)        { ny++; shifty--; }
      if (shiftz    > 0)        { nz++; shiftz--; }
      if (shiftx+nx < dims[0])  { nx++;           }
      if (shifty+ny < dims[1])  { ny++;           }
      if (shiftz+nz < dims[2])  { nz++;           }

      nodes        .resize((nx+1)*(ny+1)*(nz+1));
      local2Global_.resize((nx+1)*(ny+1)*(nz+1));
      elements     .resize(nx*ny*nz);
      
#define NODE(i,j,k) ((k)*(ny+1)*(nx+1) + (j)*(nx+1) + (i))
#define CELL(i,j,k) ((k)*ny*nx         + (j)*nx     + (i))
      for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
          for (int i = 0; i <= nx; i++) {
            int ii = shiftx+i, jj = shifty+j, kk = shiftz+k;
            int nodeID = NODE(i,j,k);
            nodes[nodeID]         = Point((ii+1)*hx_, (jj+1)*hy_, (kk+1)*hz_);
            local2Global_[nodeID] = kk*nx_*ny_ + jj*nx_ + ii;
          }

      for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
          for (int i = 0; i < nx; i++) {
            std::vector<int>& element = elements[CELL(i,j,k)];
            element.resize(8);
            element[0] = NODE(i,  j,   k  );
            element[1] = NODE(i+1,j,   k  );
            element[2] = NODE(i,  j+1, k  );
            element[3] = NODE(i+1,j+1, k  );
            element[4] = NODE(i,  j,   k+1);
            element[5] = NODE(i+1,j,   k+1);
            element[6] = NODE(i,  j+1, k+1);
            element[7] = NODE(i+1,j+1, k+1);
          }
#undef NODE
#undef CELL
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildPoints(std::vector<Point>& quadPoints, std::vector<double>& quadWeights) {

      // Use nonstandard quadrature for Helmholtz from the following paper:
      // M. Ainsworth, H.A. Wajid. Optimally blended spectral-finite element scheme for wave propagation
      // and nonstandard reduced integration. SIAM J. Numer. Analysis, 48(1): 346-371, 2010.

      // Make reference domain [0,hx_] x [0,hy_] x [0,hz_].
      double qdpt = 0.816496580927726;
      quadPoints.resize(8);
      quadPoints[0] = Point( (1.0+qdpt)*hx_/2.0, (1.0+qdpt)*hy_/2.0, (1.0+qdpt)*hz_/2.0 );
      quadPoints[1] = Point( (1.0+qdpt)*hx_/2.0, (1.0-qdpt)*hy_/2.0, (1.0+qdpt)*hz_/2.0 );
      quadPoints[2] = Point( (1.0-qdpt)*hx_/2.0, (1.0+qdpt)*hy_/2.0, (1.0+qdpt)*hz_/2.0 );
      quadPoints[3] = Point( (1.0-qdpt)*hx_/2.0, (1.0-qdpt)*hy_/2.0, (1.0+qdpt)*hz_/2.0 );
      quadPoints[4] = Point( (1.0+qdpt)*hx_/2.0, (1.0+qdpt)*hy_/2.0, (1.0-qdpt)*hz_/2.0 );
      quadPoints[5] = Point( (1.0+qdpt)*hx_/2.0, (1.0-qdpt)*hy_/2.0, (1.0-qdpt)*hz_/2.0 );
      quadPoints[6] = Point( (1.0-qdpt)*hx_/2.0, (1.0+qdpt)*hy_/2.0, (1.0-qdpt)*hz_/2.0 );
      quadPoints[7] = Point( (1.0-qdpt)*hx_/2.0, (1.0-qdpt)*hy_/2.0, (1.0-qdpt)*hz_/2.0 );
      quadWeights.resize(8);
      quadWeights[0] = hx_*hy_*hz_/8.0;
      quadWeights[1] = hx_*hy_*hz_/8.0;
      quadWeights[2] = hx_*hy_*hz_/8.0;
      quadWeights[3] = hx_*hy_*hz_/8.0;
      quadWeights[4] = hx_*hy_*hz_/8.0;
      quadWeights[5] = hx_*hy_*hz_/8.0;
      quadWeights[6] = hx_*hy_*hz_/8.0;
      quadWeights[7] = hx_*hy_*hz_/8.0;

    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalBasis(Point& quadPoint, std::vector<double>& vals, std::vector<double>& dxs, std::vector<double>& dys, std::vector<double>& dzs) {

      // For the reference element, evaluate each basis function and its derivatives at a particular quadrature point
      int numBasis=(px_+1)*(py_+1)*(pz_+1);
      vals.resize(numBasis);
      dxs.resize(numBasis);
      dys.resize(numBasis);
      dzs.resize(numBasis);
      double x=quadPoint.x;
      double y=quadPoint.y;
      double z=quadPoint.z;
      double valx1=(hx_-x)/hx_;   double dx1=-1/hx_;
      double valx2=x/hx_;         double dx2=1/hx_;
      double valy1=(hy_-y)/hy_;   double dy1=-1/hy_;
      double valy2=y/hy_;         double dy2=1/hy_;
      double valz1=(hz_-z)/hz_;   double dz1=-1/hz_;
      double valz2=z/hz_;         double dz2=1/hz_;
      // linear function values, derivatives in x, y, and z directions   
      vals[0]=valx1*valy1*valz1;  dxs[0]=dx1*valy1*valz1;  dys[0]=valx1*dy1*valz1;  dzs[0]=valx1*valy1*dz1;
      vals[1]=valx2*valy1*valz1;  dxs[1]=dx2*valy1*valz1;  dys[1]=valx2*dy1*valz1;  dzs[1]=valx2*valy1*dz1;
      vals[2]=valx1*valy2*valz1;  dxs[2]=dx1*valy2*valz1;  dys[2]=valx1*dy2*valz1;  dzs[2]=valx1*valy2*dz1;
      vals[3]=valx2*valy2*valz1;  dxs[3]=dx2*valy2*valz1;  dys[3]=valx2*dy2*valz1;  dzs[3]=valx2*valy2*dz1;
      vals[4]=valx1*valy1*valz2;  dxs[4]=dx1*valy1*valz2;  dys[4]=valx1*dy1*valz2;  dzs[4]=valx1*valy1*dz2;
      vals[5]=valx2*valy1*valz2;  dxs[5]=dx2*valy1*valz2;  dys[5]=valx2*dy1*valz2;  dzs[5]=valx2*valy1*dz2;
      vals[6]=valx1*valy2*valz2;  dxs[6]=dx1*valy2*valz2;  dys[6]=valx1*dy2*valz2;  dzs[6]=valx1*valy2*dz2;
      vals[7]=valx2*valy2*valz2;  dxs[7]=dx2*valy2*valz2;  dys[7]=valx2*dy2*valz2;  dzs[7]=valx2*valy2*dz2;

    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void HelmholtzFEM3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalStretch(double& shiftx, double& shifty, double& shiftz, const std::vector<Point>& quadPoints,
												      std::vector<Scalar>& sx, std::vector<Scalar>& sy, std::vector<Scalar>& sz) {
      // For the current element domain, evaluate the PML stretching functions at the quadrature points
      sx.resize(quadPoints.size());
      sy.resize(quadPoints.size());
      sz.resize(quadPoints.size());
      for(int i=0; i<quadPoints.size(); i++) {
	double quadx=quadPoints[i].x;
	double quady=quadPoints[i].y;
	double quadz=quadPoints[i].z;
	double curx=shiftx+quadx;
	double cury=shifty+quady;
	double curz=shiftz+quadz;
	double sigx, sigy, sigz;
	if(curx<LBx_)        { sigx = delta_*pow((curx-LBx_)/PMLwidthx_,2);  }
	else if(curx>RBx_)   { sigx = delta_*pow((curx-RBx_)/PMLwidthx_,2);  }
	else                 { sigx = 0.0;                                   }
	if(cury<LBy_)        { sigy = delta_*pow((cury-LBy_)/PMLwidthy_,2);  }
	else if(cury>RBy_)   { sigy = delta_*pow((cury-RBy_)/PMLwidthy_,2);  }
	else                 { sigy = 0.0;                                   }
	if(curz<LBz_)        { sigz = delta_*pow((curz-LBz_)/PMLwidthz_,2);  }
	else if(curz>RBz_)   { sigz = delta_*pow((curz-RBz_)/PMLwidthz_,2);  }
	else                 { sigz = 0.0;                                   }
	Scalar pmlx(1.0,sigx);
	Scalar pmly(1.0,sigy);
	Scalar pmlz(1.0,sigz);
	sx[i]=pmlx;
	sy[i]=pmly;
	sz[i]=pmlz;
      }

    }


  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_HELMHOLTZ3DPROBLEM_HPP
