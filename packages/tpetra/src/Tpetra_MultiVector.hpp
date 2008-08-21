// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// FINISH: some of these arrayview objects should be something else, like Ptr

#ifndef TPETRA_MULTIVECTOR_HPP
#define TPETRA_MULTIVECTOR_HPP


#include <Teuchos_TestForException.hpp>
#include "Tpetra_MultiVectorDecl.hpp"
#include "Tpetra_MultiVectorData.hpp"

namespace Tpetra {


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, Ordinal NumVectors, bool zeroOut) 
    : DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::MultiVector")
    , MVData_()
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source) 
    : DistObject<Ordinal,Scalar>(source)
    , MVData_(source.MVData_)
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &A, Ordinal LDA, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::MultiVector")
    , MVData_()
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &arrayOfArrays, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::MultiVector")
    , MVData_()
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source, const Teuchos::ArrayView<const Ordinal> &indices)
    : DistObject<Ordinal,Scalar>(source)
    , MVData_(source.MVData_)
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source, Ordinal startIndex, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(source)
    , MVData_(source.MVData_)
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::~MultiVector()
  {
    TEST_FOR_EXCEPT(true);
  }

  template <typename Ordinal, typename Scalar> 
  void MultiVector<Ordinal,Scalar>::print(std::ostream &os) const
  {
    TEST_FOR_EXCEPT(true);
  }
  
  template <typename Ordinal, typename Scalar> 
  void MultiVector<Ordinal,Scalar>::printValues(std::ostream &os) const
  {
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  bool MultiVector<Ordinal,Scalar>::checkSizes(const DistObject<Ordinal,Scalar> &sourceObj) 
  {
    TEST_FOR_EXCEPT(true);
    return true;
  }

  template<typename Ordinal, typename Scalar>
  int MultiVector<Ordinal,Scalar>::copyAndPermute(
      const DistObject<Ordinal,Scalar> &sourceObj,
      Ordinal numSameIDs, Ordinal numPermuteIDs,
      const std::vector<Ordinal> &permuteToLIDs, const std::vector<Ordinal> &permuteFromLIDs)
  {
    TEST_FOR_EXCEPT(true);
    return 0;
  }

  template<typename Ordinal, typename Scalar>
  int MultiVector<Ordinal,Scalar>::packAndPrepare(
      const DistObject<Ordinal,Scalar> &sourceObj,
      Ordinal numExportIDs,
      const std::vector<Ordinal> &exportLIDs, std::vector<Scalar> &exports,
      Ordinal &packetSize,
      Distributor<Ordinal> &distor) 
  {
    TEST_FOR_EXCEPT(true);
    return 0;
  }

  template<typename Ordinal, typename Scalar>
  int MultiVector<Ordinal,Scalar>::unpackAndCombine(
      Ordinal numImportIDs,
      const std::vector<Ordinal> &importLIDs,
      const std::vector<Scalar> &imports,
      Distributor<Ordinal> &distor,
      CombineMode CM) 
  {
    TEST_FOR_EXCEPT(true);
    return 0;
  }
  
  /*

      // Basic constructor
      MultiVector(const VectorSpace<Ordinal, Scalar>& vectorSpace, const Ordinal NumVectors) :
        vectorSpace_(vectorSpace),
        NumVectors_(NumVectors)
      {
        Init(); 

        array_.resize(getNumVectors());

        for (Ordinal i = OrdinalZero_ ; i < NumVectors ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<Ordinal, Scalar> (vectorSpace));
        }
      }

      // Creates a deep copy of the input set of vectors.
      MultiVector(const VectorSpace<Ordinal, Scalar>& vectorSpace, 
                  std::vector<Tpetra::Vector<Ordinal, Scalar> const *> list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init(); 

        array_.resize(getNumVectors());

        for (Ordinal i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // deep copy of each of the vectors
          array_[i] = Teuchos::rcp(new Vector<Ordinal, Scalar> (*(list[i])));
        }
      }

      // Creates a shallow copy of the input set of vectors.
      MultiVector(const VectorSpace<Ordinal, Scalar>& vectorSpace, 
                  std::vector<Teuchos::RCP<Tpetra::Vector<Ordinal, Scalar> > > list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init();

        array_.resize(NumVectors_);

        for (Ordinal i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // copy RCP's from the list to this array.
          array_[i] = list[i];
        }
      }

      // Copy constructor.
      MultiVector(const MultiVector& rhs) :
        vectorSpace_(rhs.vectorSpace()),
        NumVectors_(rhs.getNumVectors())
      {
        Init();

        array_.resize(NumVectors_);

        for (Ordinal i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<Ordinal, Scalar> (*(rhs.GetVector(i))));
        }
      }

      //! Returns the global number of entries.
      Ordinal getNumGlobalEntries() const
      {
        return(vectorSpace_.getNumGlobalEntries());
      }

      //! Returns the number of entries on the calling image.
      Ordinal getNumMyEntries() const
      {
        return(vectorSpace_.getNumMyEntries());
      }

      //! Returns the number of vectors in this multivector.
      Ordinal getNumVectors() const
      {
        return(NumVectors_);
      }

      //! Returns a reference to the vector space of this multivector.
      VectorSpace<Ordinal, Scalar> const& vectorSpace () const
      {
        return(vectorSpace_);
      }

      //! Returns a reference to the i-th element of the j-th vector.
      Scalar& operator() (const Ordinal i, const Ordinal j)
      {
        return((*array_[j])[i]);
      }

      //! Returns a reference to the i-th element of the j-th vector (const version)
      Scalar const& operator() (const Ordinal i, const Ordinal j) const
      {
        return((*array_[j])[i]);
      }

      //! Sets all elements of all vectors to the given value.
      void setAllToScalar(Scalar const value)
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          array_[i]->setAllToScalar(value);
      }

      //! Sets all elements of all vectors to random value.
      void setAllToRandom()
      {
        // FIXME: sets only the real part to random
        for (int i = 0 ; i < NumVectors_ ; ++i)
        {
          //array_[i]->setAllToRandom();
          for (int j = 0 ; j < array_[0]->getNumMyEntries() ; ++j)
          {
            // FIXME (*array_[i])[j] = complex<double>(Teuchos::ScalarTraits<double>::random(), 0.0);
            (*array_[i])[j] = Teuchos::ScalarTraits<Scalar>::random();
          }
        }
      }

      //! Prints the vector to cout. FIXME
      void Print() const
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          cout << (*array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector.
      Teuchos::RCP<Tpetra::Vector<Ordinal, Scalar> > GetRCP(const int i)
      {
        return(array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector (const version).
      Teuchos::RCP<Tpetra::Vector<Ordinal, Scalar> const > GetRCP(const int i) const
      {
        return(array_[i]);
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector.
      Tpetra::Vector<Ordinal, Scalar>* GetVector(const int i)
      {
        return(array_[i].get());
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector (const version).
      const Tpetra::Vector<Ordinal, Scalar>* GetVector(const int i) const
      {
        return(array_[i].get());
      }

      void norm1(Scalar* Values) const
      {
        for (Ordinal i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<Scalar>::magnitude(array_[i]->norm1());
        }
      }

      void norm2(typename Teuchos::ScalarTraits<Scalar>::magnitudeType* Values) const
      {
        for (Ordinal i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<Scalar>::magnitude(array_[i]->norm2());
        }
      }

      void normInf(Scalar* Values) const
      {
        for (Ordinal i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<Scalar>::magnitude(array_[i]->normInf());
        }
      }

      void dotProduct(const MultiVector<Ordinal, Scalar>& A, Scalar* Values) const
      {
        for (int v = 0 ; v < getNumVectors() ; ++v)
        {
          Values[v] = GetVector(v)->dotProduct(*(A.GetVector(v)));
        }
      }

    private:

      void Init()
      {
        OrdinalZero_ = Teuchos::ScalarTraits<Ordinal>::zero();
        OrdinalOne_  = Teuchos::ScalarTraits<Ordinal>::one();

        ScalarZero_ = Teuchos::ScalarTraits<Scalar>::zero();
        ScalarOne_  = Teuchos::ScalarTraits<Scalar>::one();
      }

      std::vector<Teuchos::RCP<Vector<Ordinal, Scalar> > > array_;
      const VectorSpace<Ordinal, Scalar>& vectorSpace_;
      Ordinal NumVectors_;

      Ordinal OrdinalZero_;
      Ordinal OrdinalOne_;

      Scalar ScalarZero_;
      Scalar ScalarOne_;
*/

} // namespace Tpetra

#endif // TPETRA_MULTIVECTOR_HPP
