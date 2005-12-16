// ==============================================================
// Simple class to extend Tpetra::Vector into Tpetra::MultiVector
// This is very preliminary! Don't use unless you know what you 
// are doing...
//
// This class stores a std::vector of RefCountPtr's to allocated
// Tpetra::Vector's. It overloads the (i, j) operator so that
// one can easily access the j-th element of the i-th vector.
// ==============================================================

namespace Tpetra {

  template<class OrdinalType, class ScalarType>
  class MultiVector {

    public:

      // Basic constructor
      MultiVector(const VectorSpace<OrdinalType, ScalarType>& vectorSpace, const OrdinalType NumVectors) :
        vectorSpace_(vectorSpace),
        NumVectors_(NumVectors)
      {
        Init(); 

        array_.resize(getNumVectors());

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<OrdinalType, ScalarType> (vectorSpace));
        }
      }

      // Creates a deep copy of the input set of vectors.
      MultiVector(const VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                  std::vector<Tpetra::Vector<OrdinalType, ScalarType> const *> list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init(); 

        array_.resize(getNumVectors());

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // deep copy of each of the vectors
          array_[i] = Teuchos::rcp(new Vector<OrdinalType, ScalarType> (*(list[i])));
        }
      }

      // Creates a shallow copy of the input set of vectors.
      MultiVector(const VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                  std::vector<Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> > > list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init();

        array_.resize(NumVectors_);

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // copy RefCountPtr's from the list to this array.
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

        for (OrdinalType i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<OrdinalType, ScalarType> (*(rhs.GetVector(i))));
        }
      }

      //! Returns the global number of entries.
      OrdinalType getNumGlobalEntries() const
      {
        return(vectorSpace_.getNumGlobalEntries());
      }

      //! Returns the number of entries on the calling image.
      OrdinalType getNumMyEntries() const
      {
        return(vectorSpace_.getNumMyEntries());
      }

      //! Returns the number of vectors in this multivector.
      OrdinalType getNumVectors() const
      {
        return(NumVectors_);
      }

      //! Returns a reference to the vector space of this multivector.
      VectorSpace<OrdinalType, ScalarType> const& vectorSpace () const
      {
        return(vectorSpace_);
      }

      //! Returns a reference to the i-th element of the j-th vector.
      ScalarType& operator() (const OrdinalType i, const OrdinalType j)
      {
        return((*array_[j])[i]);
      }

      //! Returns a reference to the i-th element of the j-th vector (const version)
      ScalarType const& operator() (const OrdinalType i, const OrdinalType j) const
      {
        return((*array_[j])[i]);
      }

      //! Sets all elements of all vectors to the given value.
      void setAllToScalar(ScalarType const value)
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
            (*array_[i])[j] = Teuchos::ScalarTraits<ScalarType>::random();
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
      Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> > GetRCP(const int i)
      {
        return(array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector (const version).
      Teuchos::RefCountPtr<Tpetra::Vector<OrdinalType, ScalarType> const > GetRCP(const int i) const
      {
        return(array_[i]);
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector.
      Tpetra::Vector<OrdinalType, ScalarType>* GetVector(const int i)
      {
        return(array_[i].get());
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector (const version).
      const Tpetra::Vector<OrdinalType, ScalarType>* GetVector(const int i) const
      {
        return(array_[i].get());
      }

      void norm1(ScalarType* Values) const
      {
        for (OrdinalType i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(array_[i]->norm1());
        }
      }

      void norm2(typename Teuchos::ScalarTraits<ScalarType>::magnitudeType* Values) const
      {
        for (OrdinalType i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(array_[i]->norm2());
        }
      }

      void normInf(ScalarType* Values) const
      {
        for (OrdinalType i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(array_[i]->normInf());
        }
      }

      void dotProduct(const MultiVector<OrdinalType, ScalarType>& A, ScalarType* Values) const
      {
        for (int v = 0 ; v < getNumVectors() ; ++v)
        {
          Values[v] = GetVector(v)->dotProduct(*(A.GetVector(v)));
        }
      }

    private:

      void Init()
      {
        OrdinalZero_ = Teuchos::ScalarTraits<OrdinalType>::zero();
        OrdinalOne_  = Teuchos::ScalarTraits<OrdinalType>::one();

        ScalarZero_ = Teuchos::ScalarTraits<ScalarType>::zero();
        ScalarOne_  = Teuchos::ScalarTraits<ScalarType>::one();
      }

      std::vector<Teuchos::RefCountPtr<Vector<OrdinalType, ScalarType> > > array_;
      const VectorSpace<OrdinalType, ScalarType>& vectorSpace_;
      OrdinalType NumVectors_;

      OrdinalType OrdinalZero_;
      OrdinalType OrdinalOne_;

      ScalarType ScalarZero_;
      ScalarType ScalarOne_;

  }; // class MultiVector

} // namespace Tpetra


