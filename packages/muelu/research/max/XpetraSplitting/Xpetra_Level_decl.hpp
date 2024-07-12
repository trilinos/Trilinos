// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Xpetra_Level_decl.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */
#ifndef XPETRA_LEVEL_DECL_HPP
#define XPETRA_LEVEL_DECL_HPP

namespace Xpetra {

template <class Scalar        = MultiVector<>::scalar_type,
          class LocalOrdinal  = typename MultiVector<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class Level {
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal> vector_type;
  typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal> multivector_type;
  typedef Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;

 public:
  //! Constructors/destructors
  //@{
  // Constructor
  Level(int, int);

  // Default destructor
  virtual ~Level(){};

  //@}

  //! Public methods
  //@{
  virtual void
  apply(const multivector_type& X, multivector_type& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const {};

  virtual bool hasTransposeApply() const { return false; }

  //! Get methods
  //@{

  //! Extract the number of regions instantiated in the level
  int GetNumRegions() const;

  GlobalOrdinal GetNumRegionNodes(GlobalOrdinal region_idx) const { return regionA_[region_idx]->getCrsGraph()->getGlobalNumRows(); };

  RCP<matrix_type> GetRegionMatrix(GlobalOrdinal region_idx) const { return regionA_[region_idx]; };

  //! Extract regionToAll from the level
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > GetRegionToAll() const;

  // Extract the level ID
  GlobalOrdinal GetLevelID() const { return levelID_; }

  //! Take a region index and returns the composite index for a given mesh node
  GlobalOrdinal GetCompositeIndex(int, GlobalOrdinal) const;

  //@}

  //@}

  //! Set up methods
  //@{
  //! Set the region galerkin operators
  void SetA(Array<RCP<matrix_type> >&);

  //! Set the region prolongators
  void SetP(Array<RCP<matrix_type> >&);

  //! Set the region restrioctions
  void SetR(Array<RCP<matrix_type> >&);

  //! Set the region smoothers
  void SetSmoother(Array<RCP<multivector_type> >&);

  //! Set the regionToAll structure for the current level by using information coming from the finer level
  void SetRegionToAll(Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > >);

  //@}

  //! Check method
  //@{
  //! Control that algebraic quantities have matching dimensions
  void checkConsistency() const;
  //@}

  //! Constrcution of region smoothers
  //@{
  void ComputeRegionJacobi();
  //@}

 private:
  //! Private variables
  //@{
  GlobalOrdinal levelID_     = -1;
  GlobalOrdinal num_regions_ = -1;

  //! Region operators
  Array<RCP<matrix_type> > regionA_;

  //! Region grid transfers
  Array<RCP<matrix_type> > regionP_;
  Array<RCP<matrix_type> > regionR_;

  //! Smoother
  Array<RCP<vector_type> > regionSmoother_;

  //! Auxiliary quantities to handle regions at each coarsening level
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > level_regionToAll_;

  //@}
};

}  // namespace Xpetra

#endif
