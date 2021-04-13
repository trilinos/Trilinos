#!/bin/bash

matlab -nodisplay -r "plotresultsC0($1); quit"; stty sane
epstopdf 'control_'$1'x'$1'.eps'
xdg-open 'control_'$1'x'$1'.pdf' &
