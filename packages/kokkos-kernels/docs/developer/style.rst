Style Guide
===========

We follow google's c++ coding style. See https://google.github.io/styleguide/cppguide.html and https://github.com/kokkos/kokkos-kernels/blob/master/.clang-format for details. 

.. code-block::
    :caption: Automate coding style via a pre-commit hook

        cat kokkos-kernels/.git/hooks/pre-commit
        for FILE in $(git diff --cached --name-only | egrep '.*\.cpp$|.*\.hpp$|.*\.h$')
        do
        if [ -e $file ]; then
            clang-format-8 -i -style=file $FILE
            git add $FILE
        fi
        done
        chmod +x kokkos-kernels/.git/hooks/pre-commit

.. code-block::
    :caption: Conditionally enable or disable formatting

        // clang-format off
        cpp code here
        // clang-format on

.. code-block::
    :caption: Instal clang-format on MacOS

        brew install clang-format-8

.. code-block::
    :caption: Instal clang-format on Ubuntu

        apt install clang-format-8