Design Considerations for TriBITS
---------------------------------

Some of the basic requirements and design goals for TriBITS are outlined in
the `TriBITS Overview`_ document.

.. ToDo: Discuss design requirements for TriBITS in more detail.

As stated in `TriBITS Dependency Handling Behaviors`_, `No circular
dependencies of any kind are allowed`_.  That is, no TriBITS package (or
its tests) can declare a dependency on a `downstream`_ package, period!  To
some, this might seem over constraining but adding support for circular
dependencies to the TriBITS system would add significant complexity and
space/time overhead and is a bad idea from a basic software engineering
perspective (see the *ADP (Acyclic Dependencies Principle)* in `Software
Engineering Packaging Principles`_).  From a versioning, building, and
change-prorogation perspective, any packages involved in a circular dependency
would need to be treated as a single software engineering package anyway so
TriBITS forces development teams to glob all of this stuff together into a
single `TriBITS Package`_ when cycles in software exist.  There are
numerous wonderful ways to break circular dependencies between packages that
are proven and well established in the SE community (for example, see [`Agile
Software Development, 2003`_]).

.. ToDo: Discuss why it is a good idea to explicitly list packages instead of
.. just searching for them.  Hint: Think performance and circular
.. dependencies!
