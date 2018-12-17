.. _fill_complete:

CrsMatrix Fill Complete
#######################

When changes to the ``CrsMatrix`` structure (if allowed) or values are finished,
the matrix must declare itself "fill complete" so that post matrix fill
operations can be performed.  This is accomplished by calling the
``CrsMatrix``\'s ``fillComplete()`` method.  This is an expensive operation,
because it both rearranges local data, and communicates in order to build
reusable communication patterns for sparse matrix-vector multiply. Attempts to
amortize the cost of this operation should be made whenever possible over many
sparse matrix-vector multiplies.

``fillComplete()`` takes two arguments:

* the domain ``Map`` (the distribution of the input vector :math:`x` in a sparse
  matrix-vector multiply :math:`y = Ax`); and
* the range ``Map`` (the distribution of the output vector :math:`y` in a sparse
  matrix-vector multiply :math:`y = Ax`).

Both the domain and range ``Map``\s must be one-to-one. That is, each global
index in the ``Map`` must be uniquely owned by one `and only one` process. These
two arguments must be supplied to ``fillComplete()`` if any of the following
conditions are met:

* When the row ``Map`` is not the same as the range Map (it can't be if the row
  ``Map`` is not one to one);
* When the domain and range ``Map``\s are not equal (e.g., if the matrix is not
  square);
* When the domain or range ``Map`` as not the same as the row ``Map``.

If the domain and range ``Map``\s are the same as the row ``Map`` and the row
``Map`` is one-to-one, then you may call ``fillComplete()`` with no arguments.

Comparison with Epetra
----------------------

The most significant difference between Epetra and Tpetra sparse matrices, is
that in order to modify the entries of a ``Tpetra::CrsMatrix`` once you have
called ``fillComplete()``, you must first call ``resumeFill()``.
``Epetra_CrsMatrix`` has no corresponding "resume fill" method, and you may
modify the values of entries after ``FillComplete()`` has been called.

The reason for this difference is that Tpetra's ``fillComplete()`` has the right
to rearrange the matrix's data in ways that violate user expectations. For
example, it may give the data to a third-party library that rearranges them in
an opaque way, or even copy them into a different memory space (for example,
into GPU memory). Calling ``resumeFill()`` signals Tpetra that you want to
change either the values or the structure.
