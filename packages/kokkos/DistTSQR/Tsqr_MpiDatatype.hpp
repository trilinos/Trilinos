#ifndef __TSQR_MpiDatatype_hpp
#define __TSQR_MpiDatatype_hpp

#include <Teuchos_ConfigDefs.hpp> // HAVE_MPI

#ifdef HAVE_MPI
#  include <mpi.h>
#  include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace MPI {

    /// \brief Clone an MPI_Datatype object.
    ///
    MPI_Datatype 
    cloneRawDatatype (MPI_Datatype in, const bool needsFree);

    /// Implements a map from C++ datatype to MPI_Datatype.
    ///
    /// MpiDatatype manages type creation and freeing (via
    /// MPI_Type_free()) automatically, for types for which it's
    /// necessary.  For other types, MpiDatatype::get() just returns
    /// the predefined constant (e.g., MPI_DOUBLE or MPI_INT).  Note
    /// also that MpiDatatype< Datum > has only been defined for
    /// certain types of Datum (see MpiDatatype.cpp for which ones).
    template< class Datum >
    class MpiDatatype {
    public:
      /// Constructor, specialized by hand for typical scalar data
      /// types.
      MpiDatatype ();

      /// Copy constructor
      ///
      MpiDatatype (const MpiDatatype& rhs) 
      {
	clone (this->type_, this->needsFree_, rhs);
      }

      /// Assignment operator
      ///
      MpiDatatype& operator= (const MpiDatatype& rhs)
      {
	if (this != &rhs)
	  {
	    if (needsFree_)
	      {
		// Return value doesn't matter...
		(void) MPI_Type_free (&type_);
		needsFree_ = false;
	      }
	    clone (this->type_, this->needsFree_, rhs);
	  }
	return *this;
      }

      /// Destructor
      ///
      ~MpiDatatype () 
      { 
	if (needsFree_)
	  {
	    // Return value doesn't matter...
	    (void) MPI_Type_free (&type_);
	    needsFree_ = false;
	  }
      }

      MPI_Datatype get() const { return type_; }

    private:
      static void 
      clone (MPI_Datatype& newType, 
	     bool& needsFree,
	     const MpiDatatype& rhs)
      {
	newType = cloneRawDatatype (rhs.get(), rhs.needsFree_);
	needsFree = rhs.needsFree_;
      }

      /// The actual MPI_Datatype object corresponding to Datum.
      ///
      MPI_Datatype type_;

      /// Whether or not MPI_Type_free() needs to be called on type_, in
      /// ~MpiDatatype().
      bool needsFree_;
    };

  } // namespace MPI
} // namespace TSQR


#endif // HAVE_MPI
#endif // __TSQR_MpiDatatype_hpp
