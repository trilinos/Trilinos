# Doxygen

`MiniTensor`'s documentation is built in the build tree:

1. Run your `cmake` configure line or script.
1. In the <buildDirectory>, type `make doc_minitensor`.

The documentation will be placed in `<buildDirectory>/packages/minitensor/doc/html`.

## Publishing to the Trilinos website

The Trilinos website serves package documentation from the
[trilinos.github.io](https://github.com/trilinos/trilinos.github.io)
repository: the generated HTML for a package `<pkg>` lives under
`docs/<pkg>/` there and is served at
`https://trilinos.github.io/docs/<pkg>/index.html`.

To publish or refresh MiniTensor's pages:

1. Build the documentation as above (`make doc_minitensor`).
1. Clone `trilinos/trilinos.github.io` and replace the contents of
   `docs/minitensor/` with the contents of
   `<buildDirectory>/packages/minitensor/doc/html/`.
1. Open a pull request against `trilinos/trilinos.github.io`.
