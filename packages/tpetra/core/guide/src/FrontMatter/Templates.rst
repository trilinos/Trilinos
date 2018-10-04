.. _templated_types_in_tpetra:

Templated Types in Tpetra
#########################

All classes in Tpetra are templated.  Templates allow:

* specifying types of ordinals (matrix/vector indices) and scalars;
* choosing types as needed for their application; and
* choosing types that increase functionality.

For example, 64-bit ordinals allow for problem sizes to break the 2 billion
element barrier present in Epetra, whereas complex scalar types allow the native
description and solution of complex-valued problems.

Standard Tpetra Templated Types
===============================

Most of the classes in Tpetra are templated according to the data types which
constitute the class. These are the following:

* ``Scalar``: A ``Scalar`` is the type of values in the sparse
  matrix or dense vector. This is the type most likely to be changed by many
  users. The most common use cases are :cpp:``float``, ``double``,
  ``std::complex<float>`` and ``std::complex<double>``.
  However, many other data types can be used, as long as they have
  specializations for ``Teuchos::ScalarTraits`` and
  ``Teuchos::SerializationTraits``, and support the necessary arithmetic
  operations, such as addition, subtraction, division and multiplication.

* ``LocalOrdinal``: A ``LocalOrdinal`` is used to store indices representing
  local IDs. The standard use case, as well as the default for most classes, is
  ``int``. Any signed built-in C++ integer type may be used. The reason why local
  and global ordinals may have different types is for efficiency. If the
  application allows it, using smaller local ordinals requires less storage and
  may improve performance of computational kernels such as sparse matrix-vector
  multiply.  The ``LocalOrdinal`` is ``typedef``'d in most Tpetra objects as
  ``local_ordinal_type``.

* ``GlobalOrdinal``: A ``GlobalOrdinal`` is used to store global indices and to
  describe global properties of a distributed object (e.g., global number of
  entries in a sparse matrix, or global number of rows in a vector.) The
  ``GlobalOrdinal``, therefore, dictates the maximum size of a distributed
  object.  Like the ``LocalOrdinal``, it should be signed built-in C++ integer
  type and its default type is ``int``. If a problem has more than 2 billion
  entries, a 64-bit integer type will be needed (such as ``long long`` or
  ``int64_t``) for ``GlobalOrdinal``. Note, however, that if you have enough
  processes so that no one process stores more than 2 billion entries locally,
  then a 32-bit integer type (such as ``int`` or ``int32_t``) can still be used
  for ``LocalOrdinal``.  ``GlobalOrdinal`` is ``typedef``'d in most Tpetra
  objects as ``global_ordinal_type``.

* ``Node``: Computational classes in Tpetra will also be templated on a ``Node``
  type. This specifies the node-level parallel programming model for Tpetra
  objects to use.

Guidance on Choice of Ordinals
==============================

It is usually more efficient to use the shortest integer type possible for both
local and global indices. "Shortest" means fewest number of bits. Fewer bits
mean you use less memory and thus can solve bigger problems or use
higher-quality preconditioners that solve problems in fewer iterations. Shorter
local indices can also mean better performance for local sparse matrix kernels,
such as sparse matrix-vector multiply, sparse triangular solve, and smoothing
(for algebraic multigrid).

Tpetra differs from Epetra in that you, the user, get to decide the types of
local and global indices. In Epetra, local and global indices both used to have
type int. With the latest Trilinos release, Epetra may use either 32-bit or
64-bit integers for global indices, and always uses 32-bit integers (int) for
local indices. Tpetra lets you decide the types of each.
