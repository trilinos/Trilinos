// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CONTAINERFACTORY_DECL_H
#define IFPACK2_CONTAINERFACTORY_DECL_H

#include "Ifpack2_Container.hpp"
#include "Ifpack2_Partitioner.hpp"
#ifdef HAVE_IFPACK2_AMESOS2
#  include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include <string>
#include <map>

/// \file Ifpack2_ContainerFactory_decl.hpp
/// \brief Ifpack2::ContainerFactory class declaration

namespace Ifpack2 {
namespace Details {

//The default container type names:
//
//Dense
//SparseILUT
//SparseAmesos, alias SparseAmesos2
//TriDi
//Banded

template<typename MatrixType>
struct ContainerFactoryEntryBase
{
  virtual Teuchos::RCP<Ifpack2::Container<MatrixType>> build(
      const Teuchos::RCP<const MatrixType>& A,
      const Teuchos::Array<Teuchos::Array<typename MatrixType::local_ordinal_type>>& partitions,
      const Teuchos::RCP<const Tpetra::Import<
        typename MatrixType::local_ordinal_type,
        typename MatrixType::global_ordinal_type,
        typename MatrixType::node_type>> importer,
      bool pointIndexed) = 0;
  virtual ~ContainerFactoryEntryBase() {}
};

template<typename MatrixType, typename ContainerType>
struct ContainerFactoryEntry : public ContainerFactoryEntryBase<MatrixType>
{
  Teuchos::RCP<Ifpack2::Container<MatrixType>> build(
      const Teuchos::RCP<const MatrixType>& A,
      const Teuchos::Array<Teuchos::Array<typename MatrixType::local_ordinal_type>>& partitions,
      const Teuchos::RCP<const Tpetra::Import<
        typename MatrixType::local_ordinal_type,
        typename MatrixType::global_ordinal_type,
        typename MatrixType::node_type>> importer,
      bool pointIndexed)
  {
    return Teuchos::rcp(new ContainerType(A, partitions, importer, pointIndexed));
  }
  ~ContainerFactoryEntry() {}
};

} // namespace Details

/// \class ContainerFactory
/// \brief A static "factory" that provides a way to register
///   and construct arbitrary Ifpack2::Container subclasses using
///   string keys.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.

template<typename MatrixType>
struct ContainerFactory
{
  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;
  //! The local_ordinal_type from the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  //! The global_ordinal_type from the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  //! The node_type from the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Tpetra::RowMatrix specialization (superclass of MatrixType)
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;
  //! Tpetra::Importer specialization for use with \c MatrixType and compatible MultiVectors.
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
  typedef Ifpack2::Container<MatrixType> BaseContainer;
  //@}

  static_assert (std::is_same<typename std::decay<MatrixType>::type, row_matrix_type>::value,
                 "MatrixType must be a Tpetra::RowMatrix specialization.");

  // \name Functions
  //@{
  //! Registers a specialization of Ifpack2::Container by binding a key (string) to it.
  /*!
    \tparam ContainerType The Container specialization to register.
    \param containerType The key to pair with ContainerType. After registering, the key can be used to construct a ContainerType.
    */
  template<typename ContainerType>
  static void registerContainer(std::string containerType);

  //! Build a specialization of Ifpack2::Container given a key that has been registered.
  /*!
    \param containerType The key for looking up the Container specialization. If this key hasn't been registered, an exception is thrown.
    \param A The problem matrix.
    \param partitions The rows that correspond to each block. The outer list contains blocks, and the inner list contains rows. In BlockRelaxation, this is retrieved from a Partitioner.
    \param importer The importer that is used to import off-process rows (used by overlapping BlockRelaxation).
    \param pointIndexed If A is a BlockCrsMatrix, whether partitions contains the indices of individual DOFs instead of nodes/blocks.
    */
  static Teuchos::RCP<BaseContainer> build(std::string containerType, const Teuchos::RCP<const MatrixType>& A,
      const Teuchos::Array<Teuchos::Array<local_ordinal_type>>& partitions, const Teuchos::RCP<const import_type> importer, bool pointIndexed);

  //! Registers a specialization of Ifpack2::Container by binding a key (string) to it.
  /*!
    \param containerType The key to deregister. If it wasn't registered before, the call has no effect.
    */
  static void deregisterContainer(std::string containerType);
  //@}

  private:
  static std::map<std::string, Teuchos::RCP<Details::ContainerFactoryEntryBase<MatrixType>>> table;
  static bool registeredDefaults;     //this will initially be false
  static void registerDefaults();
};

} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CONTAINERFACTORY_H
