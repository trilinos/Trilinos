#ifndef __TSQR_TsqrApplyType_hpp
#define __TSQR_TsqrApplyType_hpp

#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class ApplyType
  /// \brief NoTranspose, Transpose, or ConjugateTranspose
  ///
  /// An ApplyType instance represents one of three ways one can apply
  /// an implicitly stored Q factor from a QR factorization to the
  /// left of a matrix C: either as Q (NoTranspose), as Q^T
  /// (Transpose), or as Q^H (ConjugateTranspose).  Transpose and
  /// ConjugateTranspose only mean different things in complex
  /// arithmetic.  This class is a kind of "checked enum" that only
  /// allows these three values.  It knows how to go from a length-one
  /// string to the appropriate ApplyType object: "N" -> NoTranspose,
  /// "T" -> Transpose, and "H" -> ConjugateTranspose.
  class ApplyType {
  public:
    /// \brief Constructor
    ///
    /// \param op [in] One of "N", "T", or "H".  Only the first
    ///   character of op is read.
    ApplyType (const std::string& op) :
      type_ (decide_apply_type (op))
    {}

    /// \brief Whether this corresponds to (Q^T or Q^H)
    ///
    /// If op corresponds to applying Q, return false, else if op
    /// corresponds to applying Q^T or Q^H, return true.
    ///
    /// \note We lump Q^T and Q^H together because they both involve
    /// applying the Q factor pieces in the same order as they were
    /// computed in factor(), whereas applying Q involves applying
    /// those pieces in the reverse order of their computation in
    /// factor().
    bool transposed () const { 
      return type_ != NoTranspose_;
    }

    ApplyType (const ApplyType& rhs) :
      type_ (rhs.type_)
    {}

    ApplyType& operator= (const ApplyType& rhs) {
      type_ = rhs.type_;
    }

    bool operator== (const ApplyType& rhs) const {
      return type_ == rhs.type_;
    }

    bool operator!= (const ApplyType& rhs) const {
      return !(type_ == rhs.type_);
    }

    static const ApplyType NoTranspose;
    static const ApplyType Transpose;
    static const ApplyType ConjugateTranspose;

  private:
    /// Enum representing the underlying three possibilities for an
    /// ApplyType: Not transposed, Transposed, or Conjugate
    /// transposed.
    enum ApplyType_ { NoTranspose_, Transpose_, ConjugateTranspose_ };
    ///
    /// The state of this ApplyType
    ///
    ApplyType_ type_;
    ///
    /// Return true if op[0] == 'T' or 'H', false otherwise
    ///
    bool 
    decide_transposed (const std::string& op) const;

    ApplyType_
    decide_apply_type (const std::string& op) const;
  };

} // namespace TSQR

#endif // __TSQR_TsqrApplyType_hpp
