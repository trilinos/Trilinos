/*
 * Xpetra_MapExtractor.hpp
 *
 *  Created on: 08.08.2011
 *      Author: tobias
 */

#ifndef XPETRA_MAPEXTRACTOR_HPP_
#define XPETRA_MAPEXTRACTOR_HPP_

#include <map>

#include <iostream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraImport.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#endif

#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_TpetraImport.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#endif

#include <Xpetra_ImportFactory.hpp>


namespace Xpetra
{
	template <class Scalar, class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
	class MapExtractor : public Teuchos::Describable
	{
		typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> MapClass;
		typedef Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> VectorClass;
		typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVectorClass;
		typedef Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> ImportClass;

	public:

		/** \name Extract subblocks from full map */
		//@{
		virtual void ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const = 0;
    virtual void ExtractVector(Teuchos::RCP<      VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const = 0;
    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block) const = 0;
    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<      VectorClass>& full, size_t block) const = 0;

    virtual void ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const = 0;
    virtual void ExtractVector(Teuchos::RCP<      MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const = 0;
    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block) const = 0;
    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<      MultiVectorClass>& full, size_t block) const = 0;
    //@}

		//virtual Teuchos::RCP<Xpetra::Vector<LocalOrdinal,GlobalOrdinal,Node> >ExtractVector(Teuchos::RCP<const Xpetra::Vector<LocalOrdinal,GlobalOrdinal,Node> > full, int block) const {};

    /** \name Insert subblocks into full map */
    //@{
		virtual void InsertVector(Teuchos::RCP<const VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const = 0;
		virtual void InsertVector(Teuchos::RCP<      VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const = 0;

    virtual void InsertVector(Teuchos::RCP<const MultiVectorClass>& partial, size_t block, Teuchos::RCP<MultiVectorClass>& full) const = 0;
    virtual void InsertVector(Teuchos::RCP<      MultiVectorClass>& partial, size_t block, Teuchos::RCP<MultiVectorClass>& full) const = 0;

		//@}

		virtual Teuchos::RCP<VectorClass> getVector(size_t i) const = 0;
		virtual Teuchos::RCP<MultiVectorClass> getVector(size_t i, size_t numvec) const = 0;

    /** \name Maps */
    //@{

    /// number of partial maps
    virtual size_t NumMaps() const = 0;

    /// get the map
    virtual const Teuchos::RCP<const MapClass> getMap(size_t i) const = 0;

    /// the full map
    virtual const Teuchos::RCP<const MapClass> getFullMap() const = 0;

    //@}

	};
}


#endif /* XPETRA_MAPEXTRACTOR_HPP_ */
