#ifndef __TSQR_Random_NormalGenerator_hpp
#define __TSQR_Random_NormalGenerator_hpp

#include <Tsqr_Lapack.hpp>
#include <algorithm>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Random {

    template< class Ordinal, class Scalar >
    class NormalGenerator {
    private:
      static const int defaultBufferLength = 100;

    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;

      /// (Pseudo)random norma(0,1) number generator, using LAPACK's
      /// _LARNV routine, wrapped in a generator interface.
      ///
      /// \param buffer_length [in] How many entries we keep buffered at
      ///   one time.  If you know how many outputs you want, set this
      ///   accordingly, so that all the expense of generation happens
      ///   at construction.
      ///
      /// \param iseed [in] Array of four integers, representing the
      ///   seed.  See documentation of _LARNV.  In particular, the
      ///   array elements must be in [0,4095], and the last element
      ///   (iseed[3]) must be odd.
      NormalGenerator (const int iseed[4], 
		       const int buffer_length = defaultBufferLength) :
	iseed_ (4),
	buffer_ (buffer_length),
	buffer_length_ (buffer_length),
	cur_pos_ (0)
      {
	std::copy (iseed, iseed+4, iseed_.begin());
	fill_buffer ();
      }


      /// (Pseudo)random normal(0,1) number generator, using LAPACK's
      /// _LARNV routine, wrapped in a generator interface.  The
      /// four-integer seed is set to [0, 0, 0, 1], which is a valid
      /// seed and which ensures a reproducible sequence.
      ///
      /// \param buffer_length [in] How many entries we keep buffered at
      ///   one time.  If you know how many outputs you want, set this
      ///   accordingly, so that all the expense of generation happens
      ///   at construction.
      NormalGenerator (const int buffer_length = defaultBufferLength) :
	iseed_ (4),
	buffer_ (buffer_length),
	buffer_length_ (buffer_length),
	cur_pos_ (0)
      {
	iseed_[0] = 0;
	iseed_[1] = 0;
	iseed_[2] = 0;
	iseed_[3] = 1;
	fill_buffer ();
      }
      
      /// Get the next value from the buffer, generating new values if
      /// necessary.  Depending on the buffer length, the generation
      /// phase may take a while.
      Scalar operator() () { return next(); }

      /// Get the current seed.  This ca be used to restart the
      /// generator, but only if you account for the buffered values.
      void 
      getSeed (int iseed[4]) const
      {
	std::copy (iseed_.begin(), iseed_.end(), iseed);
      }

    private:
      std::vector< int > iseed_;
      std::vector< Scalar > buffer_;
      int buffer_length_, cur_pos_;
      LAPACK<int, Scalar > lapack_;

      void
      fill_buffer () 
      {
	enum { uniform_0_1 = 1, 
	       uniform_m1_1 = 2, 
	       normal_0_1 = 3 } distribution_types;
	lapack_.LARNV (normal_0_1, &iseed_[0], buffer_length_, &buffer_[0]);
      }

      Scalar 
      next () 
      { 
	// Greater-than impossible, but we check for robustness' sake.
	if (cur_pos_ >= buffer_length_) 
	  {
	    fill_buffer ();
	    cur_pos_ = 0;
	  }
	return buffer_[cur_pos_++];
      }
    };
  } // namespace Random
} // namespace TSQR


#endif // __TSQR_Random_NormalGenerator_hpp
