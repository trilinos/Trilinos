#!/bin/bash

matlab -nodisplay -r "plotresultsC0($1,$2); quit"; stty sane
epstopdf 'control_'$1'x'$2'.eps'
xdg-open 'control_'$1'x'$2'.pdf' &
