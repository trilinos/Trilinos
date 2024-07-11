// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_TPETRA_1DFEM_FUNCTORS_HPP
#define NOX_TPETRA_1DFEM_FUNCTORS_HPP

#include "Kokkos_Core.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"

template<class TpetraVectorType>
struct MeshFillFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type Scalar;
  typedef typename TpetraVectorType::local_ordinal_type LO;
  typedef typename TpetraVectorType::global_ordinal_type GO;
  typedef typename TpetraVectorType::dual_view_type::t_dev ViewType;

  ViewType coordsView_;
  const Scalar zMin_;
  const Scalar dz_;
  const GO minGID_;

  MeshFillFunctor(TpetraVectorType& coords,
                  const Scalar& zMin,
                  const Scalar& dz,
                  const GO& minGID) :
    coordsView_(coords.getLocalViewDevice(Tpetra::Access::ReadWrite)),
    zMin_(zMin),
    dz_(dz),
    minGID_(minGID)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO node) const
  {
    coordsView_(node,0) = zMin_ + static_cast<Scalar>(minGID_ + node) * dz_;
  }

};

//==================================================================

template <class size_type, class LO>
struct RowCountsFunctor
{
  const Kokkos::View<size_type*> counts_;
  const LO numMyNodes_;
  const int numProcs_;
  const int myRank_;

  RowCountsFunctor(const Kokkos::View<size_type*>& counts,
                   const LO& numMyNodes,
                   const int& numProcs,
                   const int& myRank) :
    counts_(counts),
    numMyNodes_(numMyNodes),
    numProcs_(numProcs),
    myRank_(myRank)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO localRow, std::size_t& curNumLocalEntries) const
  {
    // Add a diagonal matrix entry
    Kokkos::atomic_increment(&counts_(localRow));
    ++curNumLocalEntries;
    // Contribute a matrix entry to the previous row
    if (localRow > 0) {
      Kokkos::atomic_increment(&counts_(localRow-1));
      ++curNumLocalEntries;
    }
    // Contribute a matrix entry to the next row
    if (localRow < numMyNodes_-1) {
      Kokkos::atomic_increment(&counts_(localRow+1));
      ++curNumLocalEntries;
    }
    // MPI process to the left sends us an entry
    if ((myRank_ > 0) && (localRow == 0)) {
      Kokkos::atomic_increment(&counts_(localRow));
      ++curNumLocalEntries;
    }
    // MPI process to the right sends us an entry
    if ((myRank_ < numProcs_-1) && (localRow == numMyNodes_-1)) {
      Kokkos::atomic_increment(&counts_(localRow));
      ++curNumLocalEntries;
    }
  }

};

//==================================================================

template <class size_type, class LO>
struct RowOffsetsFunctor
{
  const Kokkos::View<size_type*> offsets_;
  const Kokkos::View<size_type*> counts_;
  const LO numMyNodes_;

  RowOffsetsFunctor(const Kokkos::View<size_type*>& offsets,
                    const Kokkos::View<size_type*>& counts,
                    const LO& numMyNodes) :
    offsets_(offsets),
    counts_(counts),
    numMyNodes_(numMyNodes)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO localRow, size_type& update, const bool final) const
  {
    if (final) {
      offsets_(localRow) = update;
    }
    if (localRow < numMyNodes_) {
      update += counts_(localRow);
    }
  }

};

//==================================================================

template <class size_type, class LO>
struct ColumnIndexCompFunctor
{
  const Kokkos::View<LO*> indices_;
  const Kokkos::View<size_type*> offsets_;
  const Kokkos::View<size_type*> counts_;
  const LO numMyNodes_;
  const int numProcs_;
  const int myRank_;

  ColumnIndexCompFunctor(const Kokkos::View<LO*> indices,
                         const Kokkos::View<size_type*>& offsets,
                         const Kokkos::View<size_type*>& counts,
                         const LO& numMyNodes,
                         const int& numProcs,
                         const int& myRank) :
    indices_(indices),
    offsets_(offsets),
    counts_(counts),
    numMyNodes_(numMyNodes),
    numProcs_(numProcs),
    myRank_(myRank)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO localRow) const
  {
    // Get local column index for diagonal entry
    LO localCol = (myRank_ > 0) ? (localRow + 1) : localRow;

    // Add a diagonal matrix entry
    size_type negOne = static_cast<size_type>(-1);
    {
      // Decrementing row counts back to 0
      size_type count = Kokkos::atomic_fetch_add(&counts_(localRow), negOne) + negOne;
      indices_(offsets_(localRow) + count) = localCol;
    }
    // Contribute a matrix entry to the previous row
    if (localRow > 0) {
      size_type count = Kokkos::atomic_fetch_add(&counts_(localRow-1), negOne) + negOne;
      indices_(offsets_(localRow-1) + count) = localCol;
    }
    // Contribute a matrix entry to the next row
    if (localRow < numMyNodes_-1) {
      size_type count = Kokkos::atomic_fetch_add(&counts_(localRow+1), negOne) + negOne;
      indices_(offsets_(localRow+1) + count) = localCol;
    }
    // MPI process to the left sends us an entry
    if ((myRank_ > 0) && (localRow == 0)) {
      size_type count = Kokkos::atomic_fetch_add(&counts_(localRow), negOne) + negOne;
      indices_(offsets_(localRow) + count) = localCol - 1;
    }
    // MPI process to the right sends us an entry
    if ((myRank_ < numProcs_-1) && (localRow == numMyNodes_-1)) {
      size_type count = Kokkos::atomic_fetch_add(&counts_(localRow), negOne) + negOne;
      indices_(offsets_(localRow) + count) = localCol + 1;
    }

  }

};

//==================================================================

// Finite Element Basis Object
template <class Scalar, class LO>
class Basis {

 public:
  // Constructor
  KOKKOS_INLINE_FUNCTION
  Basis():
    uu(0.0),
    zz(0.0),
    duu(0.0),
    eta(0.0),
    wt(0.0),
    dz(0.0),
    uuold(0.0),
    duuold(0.0)
  {}

  // Destructor
  KOKKOS_INLINE_FUNCTION
  ~Basis()
  {}

  // Calculates the values of z and u at the specified Gauss point
  KOKKOS_INLINE_FUNCTION
  void computeBasis(LO gp, Scalar* z, Scalar* u, Scalar* uold = 0) {
    if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
    if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

    // Calculate basis function and derivatives at Gauss point
    phi[0]=(1.0-eta)/2.0;
    phi[1]=(1.0+eta)/2.0;
    dphide[0]=-0.5;
    dphide[1]=0.5;

    // Caculate function and derivative approximations at GP.
    dz=0.5*(z[1]-z[0]);
    zz=0.0;
    uu=0.0;
    duu=0.0;
    uuold=0.0;
    duuold=0.0;
    for (LO i=0; i < 2; i++) {
      zz += z[i] * phi[i];
      uu += u[i] * phi[i];
      duu += u[i] * dphide[i];
      if (uold) {
        uuold += uold[i] * phi[i];
        duuold += uold[i] * dphide[i];
      }
    }
  }

 public:
  // Variables that are calculated at the Gauss point
  Scalar phi[2];
  Scalar dphide[2];
  Scalar uu;
  Scalar zz;
  Scalar duu;
  Scalar eta;
  Scalar wt;
  Scalar dz;
  // These are only needed for transient
  Scalar uuold;
  Scalar duuold;
};

//==================================================================

template <class TpetraVectorType>
struct ResidualEvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type scalar_type;
  typedef typename TpetraVectorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraVectorType::map_type map_type;
  typedef typename map_type::local_map_type local_map_type;
  typedef typename TpetraVectorType::dual_view_type dual_view_type;
  typedef typename dual_view_type::t_dev view_type;
  typedef typename view_type::const_type const_view_type;

  const view_type f_view_;
  const const_view_type x_view_;
  const const_view_type u_view_;
  const local_map_type row_map_;
  const local_map_type col_map_;
  const int myRank_;
  const scalar_type k_;
  const scalar_type p4_;

  ResidualEvaluatorFunctor(TpetraVectorType& f,
                           const TpetraVectorType& x,
                           const TpetraVectorType& u,
                           const int& myRank,
                           const typename TpetraVectorType::impl_scalar_type& k,
                           const typename TpetraVectorType::impl_scalar_type& p4) :
    f_view_(f.getLocalViewDevice(Tpetra::Access::ReadWrite)),
    x_view_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    u_view_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    row_map_(f.getMap()->getLocalMap()),
    col_map_(u.getMap()->getLocalMap()),
    myRank_(myRank),
    k_(k),
    p4_(p4)
  {}

  // Adds the contribution from element elt to the residual vector
  KOKKOS_INLINE_FUNCTION
  void operator() (const local_ordinal_type elt) const
  {
    // Get the solution and coordinates at the nodes
    scalar_type xx[2];
    xx[0] = x_view_(elt,0);
    xx[1] = x_view_(elt+1,0);

    scalar_type uu[2];
    uu[0] = u_view_(elt,0);
    uu[1] = u_view_(elt+1,0);

    Basis<scalar_type, local_ordinal_type> basis;

    // Loop Over Gauss Points
    for(local_ordinal_type gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu);

      // Loop over nodes in element
      for (local_ordinal_type i = 0; i < 2; ++i) {
        local_ordinal_type localRow = row_map_.getLocalElement(col_map_.getGlobalElement(elt+i));
        const local_ordinal_type invalid = Tpetra::Details::OrdinalTraits<local_ordinal_type>::invalid();

        if (localRow != invalid) {
          scalar_type value = basis.wt * basis.dz * (basis.uu * basis.uu * basis.phi[i] * k_
            + (basis.duu * basis.dphide[i])/(basis.dz * basis.dz));
          Kokkos::atomic_add(&f_view_(localRow,0), value);
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((myRank_ == 0) && (elt == 0)) {
      scalar_type value = u_view_(0,0) - p4_;
      Kokkos::atomic_exchange(&f_view_(0,0), value);
    }

  }

};

//==================================================================

template <class TpetraVectorType, class TpetraMatrixType>
struct JacobianEvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type scalar_type;
  typedef typename TpetraVectorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraVectorType::map_type map_type;
  typedef typename map_type::local_map_type local_map_type;
  typedef typename TpetraVectorType::dual_view_type dual_view_type;
  typedef typename dual_view_type::t_dev::const_type const_view_type;
  typedef typename TpetraMatrixType::local_matrix_device_type local_matrix_type;

  const local_matrix_type J_local_;
  const const_view_type x_view_;
  const const_view_type u_view_;
  const local_map_type row_map_;
  const local_map_type col_map_;
  const int myRank_;
  const scalar_type k_;

  JacobianEvaluatorFunctor(const TpetraMatrixType& J,
                           const TpetraVectorType& x,
                           const TpetraVectorType& u,
                           const int& myRank,
                           const typename TpetraVectorType::impl_scalar_type& k) :
    J_local_(J.getLocalMatrixDevice()),
    x_view_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    u_view_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    row_map_(J.getRowMap()->getLocalMap()),
    col_map_(J.getColMap()->getLocalMap()),
    myRank_(myRank),
    k_(k)
  {}

  // Adds the contribution from element elt to the Jacobian matrix
  KOKKOS_INLINE_FUNCTION
  void operator() (const local_ordinal_type elt) const
  {
    // Get the solution and coordinates at the nodes
    scalar_type xx[2];
    xx[0]=x_view_(elt,0);
    xx[1]=x_view_(elt+1,0);

    scalar_type uu[2];
    uu[0]=u_view_(elt,0);
    uu[1]=u_view_(elt+1,0);

    Basis<scalar_type, local_ordinal_type> basis;

    // Loop Over Gauss Points
    for(local_ordinal_type gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu);

      // Loop over nodes in element
      for (local_ordinal_type i = 0; i < 2; ++i) {
        local_ordinal_type localRow = row_map_.getLocalElement(col_map_.getGlobalElement(elt+i));
        const local_ordinal_type invalid = Tpetra::Details::OrdinalTraits<local_ordinal_type>::invalid();

        if (localRow != invalid) {
          // Loop over trial functions
          for(local_ordinal_type j = 0; j < 2; ++j) {
            const local_ordinal_type localColumn = elt + j;
            scalar_type value = basis.wt * basis.dz
              * ((basis.dphide[j]*basis.dphide[i])/(basis.dz*basis.dz)
              + 2.0*basis.uu*basis.phi[j]*basis.phi[i]*k_);

            J_local_.sumIntoValues(localRow, &localColumn, 1, &value, false, true);
          }
        }
      }

      // Correct for Dirichlet BCs
      if ((myRank_ == 0) && (elt == 0)) {
        local_ordinal_type row = 0;
        local_ordinal_type column = 0;
        scalar_type value = 1.0;
        J_local_.replaceValues(row, &column, 1, &value);
        column = 1;
        value = 0.0;
        J_local_.replaceValues(row, &column, 1, &value);
      }

    }
  }

};

//==================================================================

template <class TpetraVectorType, class TpetraMatrixType>
struct PreconditionerEvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type scalar_type;
  typedef typename TpetraVectorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraVectorType::map_type map_type;
  typedef typename map_type::local_map_type local_map_type;
  typedef typename TpetraVectorType::dual_view_type dual_view_type;
  typedef typename dual_view_type::t_dev::const_type const_view_type;
  typedef typename TpetraMatrixType::local_matrix_device_type local_matrix_type;

  const local_matrix_type M_local_;
  const const_view_type x_view_;
  const const_view_type u_view_;
  const local_map_type row_map_;
  const local_map_type col_map_;
  const int myRank_;
  const scalar_type k_;

  PreconditionerEvaluatorFunctor(const TpetraMatrixType& M,
                                 const TpetraVectorType& x,
                                 const TpetraVectorType& u,
                                 const int& myRank,
                                 const typename TpetraVectorType::impl_scalar_type& k) :
    M_local_(M.getLocalMatrixDevice()),
    x_view_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    u_view_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    row_map_(M.getRowMap()->getLocalMap()),
    col_map_(M.getColMap()->getLocalMap()),
    myRank_(myRank),
    k_(k)
  {}

  // Adds the contribution from element elt to the preconditioner matrix
  KOKKOS_INLINE_FUNCTION
  void operator() (const local_ordinal_type elt) const
  {
    // Get the solution and coordinates at the nodes
    scalar_type xx[2];
    xx[0]=x_view_(elt,0);
    xx[1]=x_view_(elt+1,0);

    scalar_type uu[2];
    uu[0]=u_view_(elt,0);
    uu[1]=u_view_(elt+1,0);

    // Loop Over Gauss Points
    for(local_ordinal_type gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      Basis<scalar_type, local_ordinal_type> basis;
      basis.computeBasis(gp, xx, uu);

       // Loop over nodes in element
      for (local_ordinal_type i = 0; i < 2; ++i) {
        local_ordinal_type localRow = row_map_.getLocalElement(col_map_.getGlobalElement(elt+i));
        const local_ordinal_type invalid = Tpetra::Details::OrdinalTraits<local_ordinal_type>::invalid();

        // Loop over trial functions
        if (localRow != invalid) {
          for(local_ordinal_type j = 0; j < 2; ++j) {
            const local_ordinal_type localColumn = elt + j;
            if (row_map_.getGlobalElement(localRow) == col_map_.getGlobalElement(localColumn)) {
              scalar_type value = basis.wt * basis.dz
                * ((basis.dphide[j]*basis.dphide[i])/(basis.dz*basis.dz)
                + 2.0*basis.uu*basis.phi[j]*basis.phi[i]*k_);
              // Do we need to set force_atomic true in this?
              M_local_.sumIntoValues(localRow, &localColumn, 1, &value);
            }
          }
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((myRank_ == 0) && (elt == 0)) {
      local_ordinal_type row = 0;
      local_ordinal_type column = 0;
      scalar_type value = 1.0;
      M_local_.replaceValues(row, &column, 1, &value);
    }

  }

};

//==================================================================

template <class TpetraVectorType>
struct DfDp2EvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type scalar_type;
  typedef typename TpetraVectorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraVectorType::map_type map_type;
  typedef typename map_type::local_map_type local_map_type;
  typedef typename TpetraVectorType::dual_view_type dual_view_type;
  typedef typename dual_view_type::t_dev view_type;
  typedef typename dual_view_type::t_dev::const_type const_view_type;

  const view_type f_view_;
  const const_view_type x_view_;
  const const_view_type u_view_;
  const local_map_type row_map_;
  const local_map_type col_map_;
  const int myRank_;
  const scalar_type k_;

  DfDp2EvaluatorFunctor(TpetraVectorType& f,
                        const TpetraVectorType& x,
                        const TpetraVectorType& u,
                        const int& myRank,
                        const typename TpetraVectorType::impl_scalar_type& k) :
    f_view_(f.getLocalViewDevice(Tpetra::Access::ReadWrite)),
    x_view_(x.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    u_view_(u.getLocalViewDevice(Tpetra::Access::ReadOnly)),
    row_map_(f.getMap()->getLocalMap()),
    col_map_(u.getMap()->getLocalMap()),
    myRank_(myRank),
    k_(k)
  {}

  // Adds the contribution from element elt to the residual vector
  KOKKOS_INLINE_FUNCTION
  void operator() (const local_ordinal_type elt) const
  {
    // Get the solution and coordinates at the nodes
    scalar_type xx[2];
    xx[0] = x_view_(elt,0);
    xx[1] = x_view_(elt+1,0);

    scalar_type uu[2];
    uu[0] = u_view_(elt,0);
    uu[1] = u_view_(elt+1,0);

    Basis<scalar_type, local_ordinal_type> basis;

    // Loop Over Gauss Points
    for(local_ordinal_type gp = 0; gp < 2; ++gp) {
      // Calculate the basis function at the gauss point
      basis.computeBasis(gp, xx, uu);

      // Loop over nodes in element
      for (local_ordinal_type i = 0; i < 2; ++i) {
        local_ordinal_type localRow = row_map_.getLocalElement(col_map_.getGlobalElement(elt+i));
        const local_ordinal_type invalid = Tpetra::Details::OrdinalTraits<local_ordinal_type>::invalid();

        if (localRow != invalid) {
          scalar_type value = basis.wt * basis.dz * (basis.uu * basis.uu * basis.phi[i]);
            // + (basis.duu * basis.dphide[i])/(basis.dz * basis.dz));
          Kokkos::atomic_add(&f_view_(localRow,0), value);
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((myRank_ == 0) && (elt == 0)) {
      scalar_type value = 0.0;
      Kokkos::atomic_exchange(&f_view_(0,0), value);
    }

  }

};

#endif
