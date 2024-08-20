// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_ELASTICITY2DPROBLEM_HPP
#define GALERI_ELASTICITY2DPROBLEM_HPP

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_Problem.hpp"
#include "Galeri_MultiVectorTraits.hpp"
#include "Galeri_XpetraUtils.hpp"

namespace Galeri {

  namespace Xpetra {

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Elasticity2DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      using RealValuedMultiVector = typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector;
      Elasticity2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) {
        E  = list.get("E", Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(1e9));
        nu = list.get("nu", Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(0.25));

        nx_ = -1;
        ny_ = -1;

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

        nDim_ = 2;
        double one = 1.0;
        stretch.push_back(list.get("stretchx", one));
        stretch.push_back(list.get("stretchy", one));

        // NOTE: -1 is because galeri counts points, not elements
        dims.push_back(nx_-1);
        dims.push_back(ny_-1);

        TEUCHOS_TEST_FOR_EXCEPTION(nx_ <= 0 || ny_ <= 0, std::logic_error, "nx and ny must be positive");
        mode_ = list.get<std::string>("mode", "plane stress");
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

      void BuildMesh();
      void BuildMaterialMatrix (Teuchos::SerialDenseMatrix<LO,SC>& D);
      void BuildReferencePoints(size_t& numRefPoints, std::vector<Point>& refPoints, size_t& numGaussPoints, std::vector<Point>& gaussPoints);
    };



    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      using Teuchos::SerialDenseMatrix;
      typedef Teuchos::ScalarTraits<SC> TST;

      BuildMesh();

      const size_t numDofPerNode   = 2;
      const size_t numNodesPerElem = 4;
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
      size_t sDim = 6;
      size_t bDim = 4;
      std::vector<SerialDenseMatrix<LO,SC> > Bs(numGaussPoints);
      std::vector<SerialDenseMatrix<LO,SC> > Ss(numGaussPoints);

      for (size_t j = 0; j < numGaussPoints; j++) {
        SerialDenseMatrix<LO,SC>& S = Ss[j];
        S.shape(Teuchos::as<LO>(sDim), Teuchos::as<LO>(nDim_));
        EvalDxi (refPoints, gaussPoints[j], S[0]);
        EvalDeta(refPoints, gaussPoints[j], S[1]);

        SerialDenseMatrix<LO,SC>& B = Bs[j];
        B.shape(Teuchos::as<LO>(bDim), numDofPerElem);

        for (size_t k = 0; k < numNodesPerElem; k++) {
          B(0, numDofPerNode*k + 0) = S(k,0);
          B(1, numDofPerNode*k + 0) = S(k,1);
          B(2, numDofPerNode*k + 1) = S(k,0);
          B(3, numDofPerNode*k + 1) = S(k,1);
        }
      }

      // Construct reordering matrix (see 6.2-9 from Cook)
      SerialDenseMatrix<LO,SC> R(D->numRows(), Teuchos::as<LO>(bDim));
      R(0,0) = R(1,3) = R(2,1) = R(2,2) = 1;

      this->A_ = MatrixTraits<Map,Matrix>::Build(this->Map_, 8*numNodesPerElem);
      this->A_->setObjectLabel(this->getObjectLabel());

      SC one = TST::one(), zero = TST::zero();
      SerialDenseMatrix<LO,SC> prevKE(numDofPerElem, numDofPerElem), prevElementNodes(numNodesPerElem, Teuchos::as<LO>(nDim_));        // cache
      for (size_t i = 0; i < elements_.size(); i++) {
        // Select nodes subvector
        SerialDenseMatrix<LO,SC> elementNodes(numNodesPerElem, Teuchos::as<LO>(nDim_));
        std::vector<LO>& elemNodes = elements_[i];
        for (size_t j = 0; j < numNodesPerElem; j++) {
          elementNodes(j,0) = nodes_[elemNodes[j]].x;
          elementNodes(j,1) = nodes_[elemNodes[j]].y;
        }

        // Check if element is a translation of the previous element
        SC xMove = elementNodes(0,0) - prevElementNodes(0,0), yMove = elementNodes(0,1) - prevElementNodes(0,1);
        typename TST::magnitudeType eps = 1e-15;         // coordinate comparison criteria
        bool recompute = false;
        {
          size_t j = 0;
          for (j = 0; j < numNodesPerElem; j++)
            if (TST::magnitude(elementNodes(j,0) - (prevElementNodes(j,0) + xMove)) > eps ||
                TST::magnitude(elementNodes(j,1) - (prevElementNodes(j,1) + yMove)) > eps)
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

            SerialDenseMatrix<LO,SC> JAC(Teuchos::as<LO>(nDim_), Teuchos::as<LO>(nDim_));

            for (size_t p = 0; p < nDim_; p++)
              for (size_t q = 0; q < nDim_; q++) {
                JAC(p,q) = zero;

                for (size_t k = 0; k < numNodesPerElem; k++)
                  JAC(p,q) += S(k,p)*elementNodes(k,q);
              }

            SC detJ = JAC(0,0)*JAC(1,1) - JAC(0,1)*JAC(1,0);

            // J2 = inv([JAC zeros(2); zeros(2) JAC])
            SerialDenseMatrix<LO,SC> J2(Teuchos::as<LO>(nDim_*nDim_),Teuchos::as<LO>(nDim_*nDim_));
            J2(0,0) = J2(2,2) =  JAC(1,1) / detJ;
            J2(0,1) = J2(2,3) = -JAC(0,1) / detJ;
            J2(1,0) = J2(3,2) = -JAC(1,0) / detJ;
            J2(1,1) = J2(3,3) =  JAC(0,0) / detJ;

            SerialDenseMatrix<LO,SC> B2(J2.numRows(), B.numCols());
            B2.multiply(Teuchos::NO_TRANS,  Teuchos::NO_TRANS,    one,  J2,   B, zero);

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

                for (size_t k = 0; k < numDofPerElem; k++)
                  KE[j0][k] = KE[k][j0] = KE[j1][k] = KE[k][j1] = zero;
                KE[j0][j0] = KE[j1][j1] = one;
              }
          } else {
            // Complex case: get rid of Dirichlet DOF
            // The case is complex because if we simply reduce the size of the matrix, it would become inconsistent
            // with maps. So, instead, we modify values of the boundary cells as if we had an additional cell close
            // to the boundary. For instance, if we have a following cell
            //  D--.
            //  |  |
            //  D--D
            // we multiply all D-D connections by 2 and multiply diagonals corresponding to D by 2 or 4, depending
            // whether it is connected to 1 or 2 other D DOFs.
            for (size_t j = 0; j < numNodesPerElem; j++)
              if (dirichlet_[elemNodes[j]]) {
                LO j0 = numDofPerNode*j, j1 = j0+1;

                for (size_t k = 0; k < numNodesPerElem; k++)
                  if ((j == k) || ((j+k) & 0x1)) {
                    // Nodes j and k are connected by an edge, or j == k
                    LO k0 = numDofPerNode*k, k1 = k0+1;
                    SC f = Teuchos::as<SC>(pow(2, Teuchos::as<int>(std::min(dirichlet_[elemNodes[j]], dirichlet_[elemNodes[k]]))));

                    KE(j0,k0) *= f; KE(j0,k1) *= f;
                    KE(j1,k0) *= f; KE(j1,k1) *= f;
                }
              }
          }
        }

        // Insert KE into the global matrix
        // NOTE: KE is symmetric, therefore it does not matter that it is in the CSC format
        for (size_t j = 0; j < numDofPerElem; j++)
          if (this->Map_->isNodeGlobalElement(elemDofs[j]))
            this->A_->insertGlobalValues(elemDofs[j], elemDofs, Teuchos::ArrayView<SC>(KE[j], numDofPerElem));
      }
      this->A_->fillComplete();

      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    RCP<typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector>
    Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildCoords() {
      // FIXME: map here is an extended map, with multiple DOF per node
      // as we cannot construct a single DOF map in Problem, we repeat the coords
      this->Coords_ = MultiVectorTraits<Map,RealValuedMultiVector>::Build(this->Map_, nDim_);

      typedef typename RealValuedMultiVector::scalar_type real_type;
      typedef Teuchos::ScalarTraits<Scalar> TST;

      Teuchos::ArrayRCP<real_type> x = this->Coords_->getDataNonConst(0);
      Teuchos::ArrayRCP<real_type> y = this->Coords_->getDataNonConst(1);

      Teuchos::ArrayView<const GO> GIDs = this->Map_->getLocalElementList();

      // NOTE: coordinates vector local ordering is consistent with that of the
      // matrix map, as it is constructed by going through GIDs and translating
      // those.
      const typename TST::magnitudeType hx = TST::magnitude(stretch[0]),
                                        hy = TST::magnitude(stretch[1]);
      for (GO p = 0; p < GIDs.size(); p += 2) { // FIXME: we assume that DOF for the same node are label consequently
        GlobalOrdinal ind = GIDs[p] >> 1;
        size_t i = ind % nx_, j = ind / nx_;

        x[p] = x[p+1] = (i+1)*hx;
        y[p] = y[p+1] = (j+1)*hy;
      }

      return this->Coords_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    RCP<MultiVector> Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildNullspace() {
      const int numVectors = 3;
      this->Nullspace_ = MultiVectorTraits<Map,MultiVector>::Build(this->Map_, numVectors);

      typedef Teuchos::ScalarTraits<Scalar> TST;
      typedef typename RealValuedMultiVector::scalar_type real_type;

      if (this->Coords_ == Teuchos::null)
        BuildCoords();

      // Teuchos::ArrayView<const GO> GIDs = this->Map_->getLocalElementList();

      size_t          numDofs = this->Map_->getLocalNumElements();
      Teuchos::ArrayRCP<real_type> x = this->Coords_->getDataNonConst(0);
      Teuchos::ArrayRCP<real_type> y = this->Coords_->getDataNonConst(1);

      SC one = TST::one();

      // NOTE: nullspace local ordering is consistent with that of the matrix
      // map, as it inherits ordering from coordinates, which is consistent.

      // Translations
      Teuchos::ArrayRCP<SC> T0 = this->Nullspace_->getDataNonConst(0), T1 = this->Nullspace_->getDataNonConst(1);
      for (size_t i = 0; i < numDofs; i += nDim_) {
        T0[i]   = one;
        T1[i+1] = one;
      }

      // Calculate center
      real_type cx = this->Coords_->getVector(0)->meanValue();
      real_type cy = this->Coords_->getVector(1)->meanValue();

      // Rotations
      Teuchos::ArrayRCP<SC> R0 = this->Nullspace_->getDataNonConst(2);
      for (size_t i = 0; i < numDofs; i += nDim_) {
        // Rotate in Y-Z Plane (around Z axis): [ -y; x]
        R0[i+0] = -(y[i]-cy);
        R0[i+1] =  (x[i]-cx);
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
    void Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMesh() {
      using Teuchos::as;

      typedef Teuchos::ScalarTraits<SC> TST;
      const typename TST::magnitudeType hx = TST::magnitude(stretch[0]),
                                        hy = TST::magnitude(stretch[1]);

      GO myPID = this->Map_->getComm()->getRank();
      GO const & one=1;
      GO mx = (this->list_).get("mx", one), my = (this->list_).get("my", one);

      GO startx, starty, endx, endy;
      Utils::getSubdomainData(dims[0], mx, myPID % mx, startx, endx);
      Utils::getSubdomainData(dims[1], my, myPID / mx, starty, endy);

      LO nx = as<LO>(endx - startx), ny = as<LO>(endy - starty);

      // Expand subdomain to do overlap
      if (startx    > 0)        { nx++; startx--; }
      if (starty    > 0)        { ny++; starty--; }
      if (startx+nx < dims[0])  { nx++;        }
      if (starty+ny < dims[1])  { ny++;        }

      nodes_       .resize((nx+1)*(ny+1));
      local2Global_.resize((nx+1)*(ny+1));
      dirichlet_   .resize((nx+1)*(ny+1));
      elements_    .resize(nx*ny);

#define NODE(i,j) ((j)*(nx+1) + (i))
#define CELL(i,j) ((j)*nx     + (i))
      // NOTE: the fact that local ordering here is not consistent with that of
      // the matrix map does not matter.  The two things that matter are:
      // local2Global_ assigns to a correct GID, and nodes_ contain correct
      // coordinates
      for (LO j = 0; j <= ny; j++)
        for (LO i = 0; i <= nx; i++) {
          GO ii = startx + i, jj = starty + j;
          LO nodeID = NODE(i,j);
          nodes_       [nodeID] = Point((ii+1)*hx, (jj+1)*hy);
          local2Global_[nodeID] = jj*nx_ + ii;

          if (ii == 0   && (this->DirichletBC_ & DIR_LEFT))   dirichlet_[nodeID]++;
          if (ii == nx_ && (this->DirichletBC_ & DIR_RIGHT))  dirichlet_[nodeID]++;
          if (jj == 0   && (this->DirichletBC_ & DIR_BOTTOM)) dirichlet_[nodeID]++;
          if (jj == ny_ && (this->DirichletBC_ & DIR_TOP))    dirichlet_[nodeID]++;
        }

      for (LO j = 0; j < ny; j++)
        for (LO i = 0; i < nx; i++) {
          std::vector<LO>& element = elements_[CELL(i,j)];
          element.resize(4);
          element[0] = NODE(i,  j);
          element[1] = NODE(i+1,j);
          element[2] = NODE(i+1,j+1);
          element[3] = NODE(i,  j+1);
        }
#undef NODE
#undef CELL
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMaterialMatrix(Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& D) {
      D.shape(3,3);

      if (!strcmp(mode_.c_str(), "plane stress")) {
        typename Teuchos::ScalarTraits<SC>::magnitudeType c = E / (1 - nu*nu);
        D(0,0) = c;           D(0,1) = c*nu;
        D(1,0) = c*nu;        D(1,1) = c;
        D(2,2) = c*(1-nu)/2;

      } else if (!strcmp(mode_.c_str(), "plane strain")) {
        typename Teuchos::ScalarTraits<SC>::magnitudeType c = E / (1 + nu) / (1 - 2*nu);
        D(0,0) = c*(1-nu);    D(0,1) = c*nu;
        D(1,0) = c*nu;        D(1,1) = c*(1-nu);
        D(2,2) = c*(1-2*nu)/2;

      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown material model for 2D");
      }
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildReferencePoints(size_t& numRefPoints, std::vector<Point>& refPoints, size_t& numGaussPoints, std::vector<Point>& gaussPoints) {
      numRefPoints   = 4;
      numGaussPoints = 4;
      refPoints  .resize(numRefPoints);
      gaussPoints.resize(numGaussPoints);

      refPoints[0] = Point(-1,-1);
      refPoints[1] = Point( 1,-1);
      refPoints[2] = Point( 1, 1);
      refPoints[3] = Point(-1, 1);

      // Gauss points (reference)
      SC sq3 = 1.0/sqrt(3);
      gaussPoints[0] = Point( sq3, sq3);
      gaussPoints[1] = Point( sq3,-sq3);
      gaussPoints[2] = Point(-sq3, sq3);
      gaussPoints[3] = Point(-sq3,-sq3);
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalDxi(const std::vector<Point>& refPoints, Point& gaussPoint, SC * dxi) {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar two = one+one;
      const Scalar quarter = one/(two+two);
      for (size_t j = 0; j < refPoints.size(); j++)
        dxi[j] = refPoints[j].x * (one + refPoints[j].y*gaussPoint.y) * quarter;
      dxi[4] = -two*gaussPoint.x;
      dxi[5] = 0.0;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    void Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::EvalDeta(const std::vector<Point>& refPoints, Point& gaussPoint, SC * deta) {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar two = one+one;
      const Scalar quarter = one/(two+two);
      for (size_t j = 0; j < refPoints.size(); j++)
        deta[j] = (one + gaussPoint.x*refPoints[j].x)*refPoints[j].y * quarter;
      deta[4] = 0.0;
      deta[5] = -two*gaussPoint.y;
    }

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_ELASTICITY2DPROBLEM_HPP
