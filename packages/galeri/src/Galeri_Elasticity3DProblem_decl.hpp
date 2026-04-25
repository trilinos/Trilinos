// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_ELASTICITY3DPROBLEM_DECL_HPP
#define GALERI_ELASTICITY3DPROBLEM_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include "Galeri_Problem.hpp"

namespace Galeri {

namespace Xpetra {

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Elasticity3DProblem : public Problem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;
  Elasticity3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);

  Teuchos::RCP<Matrix> BuildMatrix();
  Teuchos::RCP<MultiVector> BuildNullspace();
  Teuchos::RCP<RealValuedMultiVector> BuildCoords();

 private:
  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  struct Point {
    SC x, y, z;

    Point() { z = Teuchos::ScalarTraits<SC>::zero(); }
    Point(SC x_, SC y_, SC z_ = Teuchos::ScalarTraits<SC>::zero())
      : x(x_)
      , y(y_)
      , z(z_) {}
  };

  GlobalOrdinal nx_, ny_, nz_;
  size_t nDim_;
  std::vector<GO> dims;
  // NOTE: nodes correspond to a local subdomain nodes. I have to construct overlapped subdomains because
  // InsertGlobalValues in Epetra does not support inserting into rows owned by other processor
  std::vector<Point> nodes_;
  std::vector<std::vector<LO> > elements_;
  std::vector<GO> local2Global_;

  std::vector<char> dirichlet_;

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType E, nu;
  std::vector<Scalar> stretch;
  std::string mode_;

  using impl_scalar_type = typename KokkosKernels::ArithTraits<SC>::val_type;
  using KAT              = KokkosKernels::ArithTraits<impl_scalar_type>;
  using Memory2D         = Kokkos::View<impl_scalar_type**, Kokkos::HostSpace>;

  void EvalD(const std::vector<Point>& refPoints, Point& gaussPoint, Memory2D S);

  void BuildMesh();
  void BuildMaterialMatrix(Memory2D& D);
  void BuildReferencePoints(size_t& numRefPoints, std::vector<Point>& refPoints, size_t& numGaussPoints, std::vector<Point>& gaussPoints);
};

}  // namespace Xpetra

}  // namespace Galeri

#endif  // GALERI_ELASTICITY3DPROBLEM_DECL_HPP
