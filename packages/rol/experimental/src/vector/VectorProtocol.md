Vector Protocol
---------------

Instead of using inheritance from an abstract vector class, free function templates
are overloaded for concrete vector types directly. To work with XROL, a type V must
satisfy certain criteria:

1) It must specialize the traits struct

    VectorElement<V>

2) It must implement a clone() method 

    template<class V>
    std::unique_ptr<V> clone( const V& x ) {
      return std::make_unique<V>( std::move( V_constructor_arguments ) );
    }

3) The vector must minimally implement *either* the elementwise functions

    template<class V, class Function, class ...Vs>
    void eval_function( V& x, const Function &f, const Vs&... vs );

    template<class R, class V>
    auto reduce( const R& r, const V& x );

  *or* the four core functions

    template<class V>
    void set( V& x, const V& y ); // x_i <- y_i

    template<class V>
    void plus( V& x, const V& y ); // x_i <- x_i + y_i

    template<class V>
    void scale( V& x, const element_t<V> alpha ); // x_i <- alpha*x_i

    template<class V>
    void fill( V& x, const element_t<V> alpha ); // x_i <- alpha
 
  If the elementwise functions are implemented, the core functions will have
  a default implementation which utilizes them as the elementwise functions 
  are more general. Alternatively, if they elementwise methods are not 
  implemented, the converse is not possible and certain features within
  XROL will not be available in a general way. For example, imposing bound
  constraints requires nonlinear operations on vectors. A generic bound constraint
  for any V which defines elementwise operations is possible, while it would
  be necessary to specialize a bound constraint for a specific V which does not.

4) There are also optional methods which may be implemented. Some of them may 
  be composed from the above functions and will have a default implementation
  which does so, for example an axpy() operation can be composed from the 
  core functions.

  Optional with default implementation in terms of core or elementwise:

    template<class V>
    magnitude_t<V> norm( const V& x );
 
    template<class V>
    void axpy( V& x, const element_t<V> alpha, const V& y ); 

  Optional methods with no default implementation:

    basis, dual, dimension, print

  There should be traits to indicate whether vector type V overloads these optional functions 
  so that attempts to use these functions leads to a compiler error if no definition exists. 

  For example, certain finite difference checks may require an implementation of basis().

