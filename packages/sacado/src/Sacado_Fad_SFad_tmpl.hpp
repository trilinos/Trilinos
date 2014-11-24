    template <typename ValueT, int Num>
    class SFad :
      public Expr< SFadExprTag<ValueT,Num > > {

    public:

      //! Base classes
      typedef Expr< SFadExprTag<ValueT,Num > > ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SFad<T,Num> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      SFad() :
        ExprType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const int sz, const ValueT & x) :
        ExprType(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      SFad(const SFad& x) :
        ExprType(static_cast<const ExprType&>(x)) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~SFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator=(const S& v) {
        ExprType::operator=(v);
        return *this;
      }

      //! Assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator=(const SFad& x) {
        ExprType::operator=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator=(const Expr<S>& x)
      {
        ExprType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator += (const S& x) {
        ExprType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator -= (const S& x) {
        ExprType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator *= (const S& x) {
        ExprType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator /= (const S& x) {
        ExprType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator += (const SFad& x) {
        ExprType::operator+=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator -= (const SFad& x) {
        ExprType::operator-=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator *= (const SFad& x) {
        ExprType::operator*=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Division-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator /= (const SFad& x) {
        ExprType::operator/=(static_cast<const ExprType&>(x));
        return *this;
      }

       //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator += (const Expr<S>& x) {
        ExprType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator -= (const Expr<S>& x) {
        ExprType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator *= (const Expr<S>& x) {
        ExprType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator /= (const Expr<S>& x) {
        ExprType::operator/=(x);
        return *this;
      }

    }; // class SFad<ValueT,Num>

    template <typename T, int Num>
    std::ostream& operator << (std::ostream& os,
                               const Expr< SFadExprTag<T,Num> >& x) {
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }
