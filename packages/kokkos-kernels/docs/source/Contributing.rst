Contributing
############


As an open source software project, the Kokkos ecosystem in general and the Kokkos Kernels library in particular, welcome and encourage contributions from any contributors. These include contribution to the source code of the library but also the build system (CMake/Spack), the testing infrastructure (workflows, container images), the documentation and our website.

Pull Requests
=============

Contribution to Kokkos Kernels should be made by pull requests against the develop branch of repository. To ease maintenance and readability of the repository we encourage contributors to create forks of the repository and feature branches in their fork before submitting changes by pull request against the main repository. Working from a fork is a requirement for leads and developers of the project. Once the pull request is opened, multiple github actions will be automatically triggered to start the integration process for the proposed changes. To merge a pull request three types requirement need to be fullfiled

#. Procedural checks need to be successful (code formatting check, Developer Certificate of Origin, dependency review, documentation).
#. Builds and tests need to be pass (CodeQL, gitlab runners, self hosted runners).
#. A code review approval from a Kokkos Kernels developer without pending request for changes from any other Kokkos Kernels developer.

For Pull Requests to be more easily reviewed, tested and merged, it is recommended to have an associated issue where the scope of the work is discussed and the changes are designed ahead of their implemenatation. Additionally, keeping the changes in a Pull Request self contained and short will make it easier for the developers to provide timely and consise reviews. Short Pull Request are also more likely to lead to well tested code!

When proposing a change, mark yourself as assignee and if possible request reviewers directly based on their expertize. When a specific reviewer is not apparent, feel free to ask for one on the kokkos-kernels slack channel. Applying labels to the Pull Request is also a good way to communicate to reviewers your intent and to make your Pull Request more easily searchable from the github API.


Testing Policy
==============

Testing is an integral part of the development of the library and is required for the healthy growth, maintenance and refactoring of the source code. Not testing the library properly exposes it to unsustainable maintenance burden and loss of confidence from the users/customers and the community.

Kokkos Kernels testing requirements are:

#. All public APIs must be tested within a unit-test, ideally if the implemtation itself is large and broken into multiple sub-kernels, each of them should be tested independently.
#. Certain algorithms require different implementations for different inputs (small vs. large matrix, backend specific implementation), all these should be exercised in unit tests.
#. When possible and appropriate, larger integration tests should be devised to stress the combination of large number of features as would happen in user applications.
#. To ensure the performance portability of the library, performance test are required for features identified as critical in common scientific and modeling applications.
#. The above requirements are applicable to all contribution from leads, developers and external collarborators.

Reviewer Behavior
=================

- Promptly provide feedback to Pull Request when a review is requested. Rapid reviews tend to keep developer engaged and leads to faster requests resolution and merge of the proposed changes.
- Request changes that are appropriate for the scope of the Pull Request, do not ask for broad changes to improve adjacent/similar code that would put an undue burden on the contributor. Instead open an issue to document possible future work that can then be prioritized by the development team.
- Mirror communication with the Pull Request author via other media on the Pull Request itself for traceability.
- If a rapid discussion can answer a question or clarify some point, prefer direct communication with the author on slack, by video chat, etc... and summarize the outcome of the discussion.
- When appropriate help the Pull Request author to raise an issue during a developers meeting to reach a path forward for resolution.
