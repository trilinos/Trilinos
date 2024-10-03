// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_ELASTICITY3DPROBLEM_HPP
#define GALERI_ELASTICITY3DPROBLEM_HPP

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_Problem.hpp"
#include "Galeri_MultiVectorTraits.hpp"
#include "Galeri_XpetraUtils.hpp"

namespace Galeri {

  namespace Xpetra {

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Elasticity3DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      using RealValuedMultiVector = typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector;
      Elasticity3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) {

        E  = list.get("E", Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(1e9));
        nu = list.get("nu", Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(0.25));

        nx_ = -1;
        ny_ = -1;
        nz_ = -1;

        if (list.isParameter("nx")) {
          if (list.isType<int>("nx"))
            nx_ = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
          else
            nx_ = list.get<GlobalOrdinal>("nx");
        }
        if (list.isParameter("ny")) {
          if (list.isType<int>("ny"))
            ny_ = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
          else
            ny_ = list.get<GlobalOrdinal>("ny");
        }
        if (list.isParameter("nz")) {
          if (list.isType<int>("nz"))
            nz_ = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
          else
            nz_ = list.get<GlobalOrdinal>("nz");
        }

        nDim_ = 3;
        double one = 1.0;
        stretch.push_back(list.get("stretchx", one));
        stretch.push_back(list.get("stretchy", one));
        stretch.push_back(list.get("stretchz", one));

        // NOTE: -1 is because galeri counts points, not elements
        dims.push_back(nx_-1);
        dims.push_back(ny_-1);
        dims.push_back(nz_-1);

        TEUCHOS_TEST_FOR_EXCEPTION(nx_ <= 0 || ny_ <= 0 || nz_ <= 0, std::logic_error, "nx, ny and nz must be positive");
      }

      Teuchos::RCP<Matrix>                BuildMatrix();
      Teuchos::RCP<MultiVector>           BuildNullspace();
      Teuchos::RCP<RealValuedMultiVector> BuildCoords();

    private:
      typedef Scalar        SC;
      typedef LocalOrdinal  LO;
      typedef GlobalOrdinal GO;

      struct Point {
        SC x, y, z;

        Point() { z = Teuchos::ScalarTraits<SC>::zero(); }
        Point(SC x_, SC y_, SC z_ = Teuchos::ScalarTraits<SC>::zero()) : x(x_), y(y_), z(z_) { }
      };

      GlobalOrdinal                  nx_, ny_, nz_;
      size_t                         nDim_;
      std::vector<GO>                dims;
      // NOTE: nodes correspond to a local subdomain nodes. I have to construct overlapped subdomains because
      // InsertGlobalValues in Epetra does not support inserting into rows owned by other processor
      std::vector<Point>             nodes_;
      std::vector<std::vector<LO> >  elements_;
      std::vector<GO>                local2Global_;

      std::vector<char>              dirichlet_;

      typename Teuchos::ScalarTraits<Scalar>::magnitudeType  E, nu;
      std::vector<Scalar>            stretch;
      std::string                    mode_;

      void EvalDxi  (const std::vector<Point>& refPoints, Point& gaussPoint, SC * dxi);
      void EvalDeta (const std::vector<Point>& refPoints, Point& gaussPoint, SC * deta);
      void EvalDzeta(const std::vector<Point>& refPoints, Point& gaussPoint, SC * dzeta);

      void BuildMesh();
      void BuildMaterialMatrix (Teuchos::SerialDenseMatrix<LO,SC>& D);
      void BuildReferencePoints(size_t& numRefPoints, std::vector<Point>& refPoints, size_t& numGaussPoints, std::vector<Point>& gaussPoints);
    };



    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      using Teuchos::SerialDenseMatrix;

      typedef Teuchos::ScalarTraits<Scalar> TST;

      BuildMesh();

      const size_t numDofPerNode   = 3;
      const size_t numNodesPerElem = 8;
      const size_t numDofPerElem   = numNodesPerElem * numDofPerNode;

      TEUCHOS_TEST_FOR_EXCEPTION(elements_[0].size() != numNodesPerElem, std::logic_error, "Incorrect number of element vertices");

      // Material constant
      SC t = 1;

      // Material matrix
      RCP<SerialDenseMatrix<LO,SC> > D(new SerialDenseMatrix<LO,SC>);
      BuildMaterialMatrix(*D);

      // Reference element, and reference Gauss points
      size_t numRefPoints, numGaussPoints;
      std::vector<Point> refPoints, gaussPoints;
      BuildReferencePoints(numRefPoints, refPoints, numGaussPoints, gaussPoints);

      // Evaluate the B matrix for the reference element
      size_t sDim = 8;
      size_t bDim = 9;
      std::vector<SerialDenseMatrix<LO,SC> > Bs(numGaussPoints);
      std::vector<SerialDenseMatrix<LO,SC> > Ss(numGaussPoints);

      for (size_t j = 0; j < numGaussPoints; j++) {
        SerialDenseMatrix<LO,SC>& S = Ss[j];
        S.shape(sDim, nDim_);
        EvalDxi  (refPoints, gaussPoints[j], S[0]);
        EvalDeta (refPoints, gaussPoints[j], S[1]);
        EvalDzeta(refPoints, gaussPoints[j], S[2]);

        SerialDenseMatrix<LO,SC>& B = Bs[j];
        B.shape(bDim, numDofPerElem);

        for (size_t k = 0; k < numNodesPerElem; k++) {
          B(0, numDofPerNode*k + 0) = S(k,0);
          B(1, numDofPerNode*k + 0) = S(k,1);
          B(2, numDofPerNode*k + 0) = S(k,2);
          B(3, numDofPerNode*k + 1) = S(k,0);
          B(4, numDofPerNode*k + 1) = S(k,1);
          B(5, numDofPerNode*k + 1) = S(k,2);
          B(6, numDofPerNode*k + 2) = S(k,0);
          B(7, numDofPerNode*k + 2) = S(k,1);
          B(8, numDofPerNode*k + 2) = S(k,2);
        }
      }

      // Construct reordering matrix (see 6.2-9 from Cook)
      SerialDenseMatrix<LO,SC> R(D->numRows(), bDim);
      R(0,0) = R(1,4) = R(2,8) = R(3,1) = R(3,3) = R(4,5) = R(4,7) = R(5,2) = R(5,6) = 1;

      this->A_ = MatrixTraits<Map,Matrix>::Build(this->Map_, numNodesPerElem * 8 * numDofPerElem);
      this->A_->setObjectLabel(this->getObjectLabel());

      SC one = Teuchos::ScalarTraits<SC>::one(), zero = Teuchos::ScalarTraits<SC>::zero();
      SerialDenseMatrix<LO,SC> prevKE(numDofPerElem, numDofPerElem), prevElementNodes(numNodesPerElem, nDim_);        // cache
      for (size_t i = 0; i < elements_.size(); i++) {
        // Select nodes subvector
        SerialDenseMatrix<LO,SC> elementNodes(numNodesPerElem, nDim_);
        std::vector<LO>& elemNodes = elements_[i];
        for (size_t j = 0; j < numNodesPerElem; j++) {
          elementNodes(j,0) = nodes_[elemNodes[j]].x;
          elementNodes(j,1) = nodes_[elemNodes[j]].y;
          elementNodes(j,2) = nodes_[elemNodes[j]].z;
        }

        // Check if element is a translation of the previous element
        SC xMove = elementNodes(0,0) - prevElementNodes(0,0), yMove = elementNodes(0,1) - prevElementNodes(0,1), zMove = elementNodes(0,2) - prevElementNodes(0,2);
        typename TST::magnitudeType eps = 1e-15;         // coordinate comparison criteria
        bool recompute = false;
        {
          size_t j = 0;
          for (j = 0; j < numNodesPerElem; j++)
            if (Teuchos::ScalarTraits<SC>::magnitude(elementNodes(j,0) - (prevElementNodes(j,0) + xMove)) > eps ||
                Teuchos::ScalarTraits<SC>::magnitude(elementNodes(j,1) - (prevElementNodes(j,1) + yMove)) > eps ||
                Teuchos::ScalarTraits<SC>::magnitude(elementNodes(j,2) - (prevElementNodes(j,2) + zMove)) > eps)
              break;
          if (j != numNodesPerElem)
            recompute = true;
        }

        SerialDenseMatrix<LO,SC> KE(numDofPerElem, numDofPerElem);
        if (recompute == false) {
          // If an element has the same form as previous element, reuse stiffness matrix
          KE = prevKE;

        } else {
          // Evaluate new stiffness matrix for the element
          SerialDenseMatrix<LO,SC> K0(D->numRows(), numDofPerElem);
          for (size_t j = 0; j < numGaussPoints; j++) {
            SerialDenseMatrix<LO,SC>& B = Bs[j];
            SerialDenseMatrix<LO,SC>& S = Ss[j];

            SerialDenseMatrix<LO,SC> JAC(nDim_, nDim_);

            for (size_t p = 0; p < nDim_; p++)
              for (size_t q = 0; q < nDim_; q++) {
                JAC(p,q) = zero;

                for (size_t k = 0; k < numNodesPerElem; k++)
                  JAC(p,q) += S(k,p)*elementNodes(k,q);
              }

            SC detJ = JAC(0,0)*JAC(1,1)*JAC(2,2) + JAC(2,0)*JAC(0,1)*JAC(1,2) + JAC(0,2)*JAC(2,1)*JAC(1,0) -
                JAC(2,0)*JAC(1,1)*JAC(0,2) - JAC(0,0)*JAC(2,1)*JAC(1,2) - JAC(2,2)*JAC(0,1)*JAC(1,0);

            // J2 = inv([JAC zeros(3) zeros(3); zeros(3) JAC zeros(3); zeros(3) zeros(3) JAC])
            SerialDenseMatrix<LO,SC> J2(nDim_*nDim_,nDim_*nDim_);
            J2(0,0) = J2(3,3) = J2(6,6) =  (JAC(2,2)*JAC(1,1)-JAC(2,1)*JAC(1,2))/detJ;
            J2(0,1) = J2(3,4) = J2(6,7) = -(JAC(2,2)*JAC(0,1)-JAC(2,1)*JAC(0,2))/detJ;
            J2(0,2) = J2(3,5) = J2(6,8) =  (JAC(1,2)*JAC(0,1)-JAC(1,1)*JAC(0,2))/detJ;
            J2(1,0) = J2(4,3) = J2(7,6) = -(JAC(2,2)*JAC(1,0)-JAC(2,0)*JAC(1,2))/detJ;
            J2(1,1) = J2(4,4) = J2(7,7) =  (JAC(2,2)*JAC(0,0)-JAC(2,0)*JAC(0,2))/detJ;
            J2(1,2) = J2(4,5) = J2(7,8) = -(JAC(1,2)*JAC(0,0)-JAC(1,0)*JAC(0,2))/detJ;
            J2(2,0) = J2(5,3) = J2(8,6) =  (JAC(2,1)*JAC(1,0)-JAC(2,0)*JAC(1,1))/detJ;
            J2(2,1) = J2(5,4) = J2(8,7) = -(JAC(2,1)*JAC(0,0)-JAC(2,0)*JAC(0,1))/detJ;
            J2(2,2) = J2(5,5) = J2(8,8) =  (JAC(1,1)*JAC(0,0)-JAC(1,0)*JAC(0,1))/detJ;

            SerialDenseMatrix<LO,SC> B2(J2.numRows(), B.numCols());
            B2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(), J2, B, zero);

            // KE = KE + t * J2B' * D * J2B * detJ
            SerialDenseMatrix<LO,SC> J2B(R.numRows(), B2.numCols());
            J2B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,    one,   R,  B2, zero);
            K0 .multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,    one,  *D, J2B, zero);
            KE .multiply(Teuchos::TRANS,    Teuchos::NO_TRANS, t*detJ, J2B,  K0, one);
          }

          // Cache the matrix and nodes
          prevKE           = KE;
          prevElementNodes = elementNodes;
        }

        Teuchos::Array<GO> elemDofs(numDofPerElem);
        for (size_t j = 0; j < numNodesPerElem; j++) { // FIXME: this may be inconsistent with the map
          elemDofs[numDofPerNode*j + 0] = local2Global_[elemNodes[j]]*numDofPerNode;
          elemDofs[numDofPerNode*j + 1] = elemDofs[numDofPerNode*j + 0] + 1;
          elemDofs[numDofPerNode*j + 2] = elemDofs[numDofPerNode*j + 0] + 2;
        }

        // Deal with Dirichlet nodes
        bool isDirichlet = false;
        for (size_t j = 0; j < numNodesPerElem; j++)
          if (dirichlet_[elemNodes[j]])
            isDirichlet = true;

        if (isDirichlet) {
          bool keepBCs = this->list_.get("keepBCs", false);
          if (keepBCs) {
            // Simple case: keep Dirichlet DOF
            // We rewrite rows and columns corresponding to Dirichlet DOF with zeros
            // The diagonal elements corresponding to Dirichlet DOF are set to 1.
            for (size_t j = 0; j < numNodesPerElem; j++)
              if (dirichlet_[elemNodes[j]]) {
                LO j0 = numDofPerNode*j+0;
                LO j1 = numDofPerNode*j+1;
                LO j2 = numDofPerNode*j+2;

                for (size_t k = 0; k < numDofPerElem; k++)
                  KE[j0][k] = KE[k][j0] = KE[j1][k] = KE[k][j1] = KE[j2][k] = KE[k][j2] = zero;
                KE[j0][j0] = KE[j1][j1] = KE[j2][j2] = one;
              }

          } else {
            // Complex case: get rid of Dirichlet DOF
            // The case is complex because if we simply reduce the size of the matrix, it would become inconsistent
            // with maps. So, instead, we modify values of the boundary cells as if we had an additional cell close
            // to the boundary.
            for (int j = 0; j < (int)numNodesPerElem; j++)
              if (dirichlet_[elemNodes[j]]) {
                LO j0 = numDofPerNode*j, j1 = j0+1, j2 = j0+2;

                // NOTE: had to make j & k int instead of size_t so that I can use subtraction without overflowing
                for (int k = 0; k < (int)numNodesPerElem; k++)
                  if ((j == k) || (std::abs(j-k) <  4 && ((j+k) & 0x1)) || (std::abs(j-k) == 4)) {
                    // Nodes j and k are connected by an edge, or j == k
                    LO k0 = numDofPerNode*k, k1 = k0+1, k2 = k0+2;
                    SC f = Teuchos::as<SC>(pow(2, Teuchos::as<int>(std::min(dirichlet_[elemNodes[j]], dirichlet_[elemNodes[k]]))));

                    KE(j0,k0) *= f; KE(j0,k1) *= f; KE(j0,k2) *= f;
                    KE(j1,k0) *= f; KE(j1,k1) *= f; KE(j1,k2) *= f;
                    KE(j2,k0) *= f; KE(j2,k1) *= f; KE(j2,k2) *= f;
                }
              }
          }
        }

        // Insert KE into the global matrix
        // NOTE: KE is symmetric, therefore it does not matter that it is in the CSC format
        for (size_t j = 0; j < numDofPerElem; j++)
          if (this->Map_->isNodeGlobalElement(elemDofs[j]))
          {
            this->A_->insertGlobalValues(elemDofs[j], elemDofs, Teuchos::ArrayView<SC>(KE[j], numDofPerElem));
          }
      }
      this->A_->fillComplete();

      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    RCP<typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector>
    Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildCoords() {
      // FIXME: map here is an extended map, with multiple DOF per node
      // as we cannot construct a single DOF map in Problem, we repeat the coords
      this->Coords_ = MultiVectorTraits<Map,RealValuedMultiVector>::Build(this->Map_, nDim_);

      typedef typename RealValuedMultiVector::scalar_type real_type;
      typedef Teuchos::ScalarTraits<Scalar> TST;

      Teuchos::ArrayRCP<real_type> x = this->Coords_->getDataNonConst(0);
      Teuchos::ArrayRCP<real_type> y = this->Coords_->getDataNonConst(1);
      Teuchos::ArrayRCP<real_type> z = this->Coords_->getDataNonConst(2);

      Teuchos::ArrayView<const GO> GIDs = this->Map_->getLocalElementList();

      // NOTE: coordinates vector local ordering is consistent with that of the
      // matrix map, as it is constructed by going through GIDs and translating
      // those.
      const typename TST::magnitudeType hx = TST::magnitude(stretch[0]),
                                        hy = TST::magnitude(stretch[1]),
                                        hz = TST::magnitude(stretch[2]);
      for (GO p = 0; p < GIDs.size(); p += 3) { // FIXME: we assume that DOF for the same node are label consequently
        GlobalOrdinal ind = GIDs[p] / 3;
        size_t i = ind % nx_, k = ind / (nx_*ny_), j = (ind - k*nx_*ny_) / nx_;

        x[p] = x[p+1] = x[p+2] = (i+1)*hx;
        y[p] = y[p+1] = y[p+2] = (j+1)*hy;
        z[p] = z[p+1] = z[p+2] = (k+1)*hz;
      }

      return this->Coords_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    RCP<MultiVector> Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildNullspace() {

      typedef Teuchos::ScalarTraits<Scalar> TST;
      typedef typename RealValuedMultiVector::scalar_type real_type;

      const int numVectors = 6;
      this->Nullspace_ = MultiVectorTraits<Map,MultiVector>::Build(this->Map_, numVectors);

      if (this->Coords_ == Teuchos::null)
        BuildCoords();

      // Teuchos::ArrayView<const GO> GIDs = this->Map_->getLocalElementList();

      size_t          numDofs = this->Map_->getLocalNumElements();
      Teuchos::ArrayRCP<real_type> x = this->Coords_->getDataNonConst(0);
      Teuchos::ArrayRCP<real_type> y = this->Coords_->getDataNonConst(1);
      Teuchos::ArrayRCP<real_type> z = this->Coords_->getDataNonConst(2);

      SC one = TST::one();

      // NOTE: nullspace local ordering is consistent with that of the matrix
      // map, as it inherits ordering from coordinates, which is consistent.

      {
	// Translations
	Teuchos::ArrayRCP<SC> T0 = this->Nullspace_->getDataNonConst(0), T1 = this->Nullspace_->getDataNonConst(1), T2 = this->Nullspace_->getDataNonConst(2);
	for (size_t i = 0; i < numDofs; i += nDim_) {
	  T0[i]   = one;
	  T1[i+1] = one;
	  T2[i+2] = one;
	}

	// Calculate center
	real_type cx = this->Coords_->getVector(0)->meanValue();
	real_type cy = this->Coords_->getVector(1)->meanValue();
	real_type cz = this->Coords_->getVector(2)->meanValue();

	// Rotations
	Teuchos::ArrayRCP<SC> R0 = this->Nullspace_->getDataNonConst(3), R1 = this->Nullspace_->getDataNonConst(4), R2 = this->Nullspace_->getDataNonConst(5);
	for (size_t i = 0; i < numDofs; i += nDim_) {
	  // Rotate in Y-Z Plane (around Z axis): [ -y; x]
	  R0[i+0] = -(y[i] - cy);
	  R0[i+1] =  (x[i] - cx);
	  
	  // Rotate in Y-Z Plane (around Z axis): [ -z; y]
	  R1[i+1] = -(z[i] - cz);
	  R1[i+2] =  (y[i] - cy);

	  // Rotate in Y-Z Plane (around Z axis): [ z; -x]
	  R2[i+0] =  (z[i] - cz);
	  R2[i+2] = -(x[i] - cx);
	}
      }

      // Equalize norms of all vectors to that of the first one
      // We do not normalize them as a vector of ones seems nice
      Teuchos::Array<typename TST::magnitudeType> norms2(numVectors);
      this->Nullspace_->norm2(norms2);
      Teuchos::Array<SC> norms2scalar(numVectors);
      for (int i = 0; i < numVectors; i++)
        norms2scalar[i] = norms2[0] / norms2[i];
      this->Nullspace_->scale(norms2scalar);

      return this->Nullspace_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMesh() {
      typedef Teuchos::ScalarTraits<Scalar> TST;
      const typename TST::magnitudeType hx = TST::magnitude(stretch[0]),
                                        hy = TST::magnitude(stretch[1]),
                                        hz = TST::magnitude(stretch[2]);

      GO myPID = this->Map_->getComm()->getRank();
      GO const & negOne=-1;
      GO mx = this->list_.get("mx", negOne), my = this->list_.get("my", negOne), mz = this->list_.get("mz", negOne);

      const GO mxy = mx*my;

      GO startx, starty, startz, endx, endy, endz;
      Utils::getSubdomainData(dims[0], mx, (myPID % mxy) % mx, startx, endx);
      Utils::getSubdomainData(dims[1], my, (myPID % mxy) / mx, starty, endy);
      Utils::getSubdomainData(dims[2], mz,  myPID / mxy,       startz, endz);

      LO nx = endx - startx, ny = endy - starty, nz = endz - startz;

      // Expand subdomain to do overlap
      if (startx    > 0)        { nx++; startx--; }
      if (starty    > 0)        { ny++; starty--; }
      if (startz    > 0)        { nz++; startz--; }
      if (startx+nx < dims[0])  { nx++;           }
      if (starty+ny < dims[1])  { ny++;           }
      if (startz+nz < dims[2])  { nz++;           }

      nodes_       .resize((nx+1)*(ny+1)*(nz+1));
      dirichlet_   .resize((nx+1)*(ny+1)*(nz+1), 0);
      local2Global_.resize((nx+1)*(ny+1)*(nz+1));
      elements_    .resize(nx*ny*nz);

#define NODE(i,j,k) ((k)*(ny+1)*(nx+1) + (j)*(nx+1) + (i))
#define CELL(i,j,k) ((k)*ny*nx         + (j)*nx     + (i))
      // NOTE: the fact that local ordering here is not consistent with that of
      // the matrix map does not matter.  The two things that matter are:
      // local2Global_ assigns to a correct GID, and nodes_ contain correct
      // coordinates
      for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
          for (int i = 0; i <= nx; i++) {
            int ii = startx+i, jj = starty+j, kk = startz+k;
            int nodeID = NODE(i,j,k);
            nodes_[nodeID]        = Point((ii+1)*hx, (jj+1)*hy, (kk+1)*hz);
            local2Global_[nodeID] = kk*nx_*ny_ + jj*nx_ + ii;

            if (ii == 0   && (this->DirichletBC_ & DIR_LEFT))   dirichlet_[nodeID]++;
            if (ii == nx_ && (this->DirichletBC_ & DIR_RIGHT))  dirichlet_[nodeID]++;
            if (jj == 0   && (this->DirichletBC_ & DIR_FRONT))  dirichlet_[nodeID]++;
            if (jj == ny_ && (this->DirichletBC_ & DIR_BACK))   dirichlet_[nodeID]++;
            if (kk == 0   && (this->DirichletBC_ & DIR_BOTTOM)) dirichlet_[nodeID]++;
            if (kk == nz_ && (this->DirichletBC_ & DIR_TOP))    dirichlet_[nodeID]++;
          }

      for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
          for (int i = 0; i < nx; i++) {
            std::vector<LO>& element = elements_[CELL(i,j,k)];
            element.resize(8);
            element[0] = NODE(i,  j,   k  );
            element[1] = NODE(i+1,j,   k  );
            element[2] = NODE(i+1,j+1, k  );
            element[3] = NODE(i,  j+1, k  );
            element[4] = NODE(i,  j,   k+1);
            element[5] = NODE(i+1,j,   k+1);
            element[6] = NODE(i+1,j+1, k+1);
            element[7] = NODE(i,  j+1, k+1);
          }
#undef NODE
#undef CELL
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMaterialMatrix(Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& D) {
      D.shape(6,6);
      typename Teuchos::ScalarTraits<SC>::magnitudeType c = E / (1 + nu) / (1 - 2*nu);
      D(0,0) = c*(1-nu);      D(0,1) = c*nu;      D(0,2) = c*nu;
      D(1,0) = c*nu;          D(1,1) = c*(1-nu);  D(1,2) = c*nu;
      D(2,0) = c*nu;          D(2,1) = c*nu;      D(2,2) = c*(1-nu);
      D(3,3) = D(4,4) = D(5,5) = c*(1-2*nu)/2;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildReferencePoints(size_t& numRefPoints, std::vector<Point>& refPoints, size_t& numGaussPoints, std::vector<Point>& gaussPoints) {
      numRefPoints   = 8;
      numGaussPoints = 8;
      refPoints  .resize(numRefPoints);
      gaussPoints.resize(numGaussPoints);

      refPoints[0] = Point(-1,-1,-1);
      refPoints[1] = Point( 1,-1,-1);
      refPoints[2] = Point( 1, 1,-1);
      refPoints[3] = Point(-1, 1,-1);
      refPoints[4] = Point(-1,-1, 1);
      refPoints[5] = Point( 1,-1, 1);
      refPoints[6] = Point( 1, 1, 1);
      refPoints[7] = Point(-1, 1, 1);

      // Gauss points (reference)
      SC sq3 = 1.0/sqrt(3);
      gaussPoints[0] = Point( sq3, sq3, sq3);
      gaussPoints[1] = Point( sq3,-sq3, sq3);
      gaussPoints[2] = Point(-sq3, sq3, sq3);
      gaussPoints[3] = Point(-sq3,-sq3, sq3);
      gaussPoints[4] = Point( sq3, sq3,-sq3);
      gaussPoints[5] = Point( sq3,-sq3,-sq3);
      gaussPoints[6] = Point(-sq3, sq3,-sq3);
      gaussPoints[7] = Point(-sq3,-sq3,-sq3);
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalDxi(const std::vector<Point>& refPoints, Point& gaussPoint, SC * dxi) {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar eight = one+one+one+one+one+one+one+one;
      const Scalar eighth = one / eight;
      for (size_t j = 0; j < refPoints.size(); j++)
        dxi[j] = refPoints[j].x * (one + refPoints[j].y*gaussPoint.y) * (one + refPoints[j].z*gaussPoint.z) * eighth;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalDeta(const std::vector<Point>& refPoints, Point& gaussPoint, SC * deta) {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar eight = one+one+one+one+one+one+one+one;
      const Scalar eighth = one / eight;
      for (size_t j = 0; j < refPoints.size(); j++)
        deta[j] = (one + refPoints[j].x*gaussPoint.x) * refPoints[j].y * (one + refPoints[j].z*gaussPoint.z) * eighth;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalDzeta(const std::vector<Point>& refPoints, Point& gaussPoint, SC * dzeta) {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar eight = one+one+one+one+one+one+one+one;
      const Scalar eighth = one / eight;
      for (size_t j = 0; j < refPoints.size(); j++)
        dzeta[j] = (one + refPoints[j].x*gaussPoint.x) * (one + refPoints[j].y*gaussPoint.y) * refPoints[j].z * eighth;
    }

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_ELASTICITY3DPROBLEM_HPP
