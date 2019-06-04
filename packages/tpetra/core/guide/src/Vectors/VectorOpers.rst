Vector Interface
################

A small subset of the |tpetra_mv|_ interface is described below.  For the full interface and more details, see the full `Tpetra Vector
documentation <tpetra_mv>`_.

Vector Template Parameters
==========================

``Tpetra::Vector`` takes four template parameters: ``Scalar``, ``LocalOrdinal``,
``GobalOrdinal``, and ``Node``.  In the following examples, default values for
the template parameters are used:

.. code-block:: c++

   typedef Tpetra::Vector<> vector_type;

The default ``Scalar``, ``LocalOrdinal``, and ``GlobalOrdinal`` types are:

.. code-block:: c++

   typedef typename vector_type::scalar_type scalar_type;
   typedef typename vector_type::local_ordinal_type local_ordinal_type;
   typedef typename vector_type::global_ordinal_type global_ordinal_type;

Vector Construction
===================

Create a Vector From a Map
--------------------------

Create a ``Vector`` from a ``Tpetra::Map``.  This version of the constructor will
fill in the vector with zeros.

.. code-block:: c++

   vector_type x(map);

where, ``map`` is a valid ``Tpetra::Map``.

With :code:`false` as the second argument, the ``Vectors``\'s data is uninitialized so
that it can be filled later without paying the cost of initially filling it with
zeros.

.. code-block:: c++

   vector_type z(map2, false);

where, ``map2`` is also a valid ``Tpetra::Map``.

Copy Constructor
----------------

Create a copy of an existing ``Vector`` using the two-argument copy constructor
with second argument ``Teuchos::Copy``.  This constructor performs a deep copy
while the one-argument copy constructor does a `shallow` copy.  ``x`` and ``y``
have the same ``Map``.

.. code-block:: c++

   vector_type y(x, Teuchos::Copy);

Vector Fill
===========

Fill with Random Entries
------------------------

Set the entries of ``x`` to (pseudo)random numbers.  Don't consider this a good
parallel pseudorandom number generator.

.. code-block:: c++

   x.randomize();

Fill with Constant Values
-------------------------

Set the entries of ``x`` to all ones

.. code-block:: c++

   x.putScalar(Teuchos::ScalarTraits<scalar_type>::one());

Note that if the conversion from ``float`` to ``scalar_type`` is defined (e.g.,
if ``scalar_type=double``), then the following is also acceptable:

.. code-block:: c++

   x.putScalar (1.0);

Vector Operations
=================

``Tpetra::Vectors``\'s interface contains some common linear algebra operations for vector-vector operations, including operations analogous to those in the BLAS 1 standard.

Vector Addition and Update
--------------------------

Set :math:`x = \beta x + \alpha z`

.. code-block:: c++

   const scalar_type alpha = 3.14159;
   const scalar_type beta = 2.71828;
   x.update (alpha, z, beta);

As long as ``x`` and ``z``\'s ``Map``\s are compatible, this is a legal
operation.  Whether it makes sense or not depends on your application.

Set :math:`y = \gamma y + \alpha x + \beta z`

.. code-block:: c++

   y.putScalar (42.0);
   const scalar_type alpha = 3.14159;
   const scalar_type beta = 2.71828;
   const scalar_type gamma = -10.0;
   y.update (alpha, x, beta, z, gamma);

Compute the 2-norm of a Vector
------------------------------

The 2-norm of a vector ``y``, :math:`\lVert y \rVert = \sqrt{\sum_{i=1}^N
y_i y_i}`,  may have a different type than ``scalar_type``.  For example, if
``scalar_type`` is complex, then the norm is real.  ``Tpetra::MultiVector`` and
``Tpetra::Vector`` give us the type of the norm.

.. code-block:: c++

  typedef typename vector_type::mag_type mag_type;
  const mag_type theNorm = y.norm2();

Scale a Vector
--------------

Scale a vector ``q`` such that :math:`q = \alpha z`

.. code-block:: c++

  q.scale(alpha, z); // q := alpha * z

.. include:: /links.txt
