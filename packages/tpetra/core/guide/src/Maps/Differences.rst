Comparison with Epetra's Maps
#############################

Tpetra's ``Map``\s look different than Epetra's because of the template
parameters, but they work similiarly. One difference is that Tpetra ``Map``\s
tend to be handled by ``Teuchos::RCP`` (reference-counted smart pointer) rather
than copied or passed by const reference. Another difference is that
``Epetra_Map`` inherits from ``Epetra_BlockMap``, whereas Tpetra's ``Map`` does
not inherit from a corresponding block map class. ``Epetra_Map`` provides only
the ``SameAs()`` predicate, whereas Tpetra's ``Map`` class distinguishes between
"compatibility" and "sameness". Finally, ``Epetra_Map``'s ``SameAs()`` is
functionally comparable to Tpetra's ``isSameAs()``.
