[![Documentation Status](https://readthedocs.org/projects/first-steps/badge/?version=latest)](https://first-steps.readthedocs.io/en/latest/?badge=latest)

# first_steps_readthedocs
First steps with readthedocs

## Tools to install

1. `pip install sphinx`
1. `pip3 install sphinx_rtd_theme`

> Note: `pip` might have to be replaced by `pip3` for `python3`.

## Setup documentation via `readthedocs`

1. `cd <path/to/project>
1. `mkdir docs`
1. `cd docs`
1. `sphinx-quickstart`

Add option
```
master_doc = 'index'
```
to `conf.py` to point the readthedocs site generator to `index.rst`, when building and deploying the site.

## Writing the documentation

Some remarks in no particular order:

- Some places require emtpy lines for proper output formatting, e.g. before a code block.
- Additional `rst` files need to be added to the `toctree` in `index.rst` to be included into the final doc site.
- Default indentation width is 3 spaces.

## Generating the documentation

- To gerenate the `html`-site, type `make html`.
- To generate a PDF via LaTeX, type `make latexpdf`.
