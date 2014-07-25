Usage:

# analyze ML prolongator smoothing using files of the form "run_*/screen_1.ml"

./arcadia.py -l MLMMtimelines -a -o screen_1 --petra=ml
Analysis is performed for screen_1.ml
                    :  i&x-mltime     eff  serialCore-mltime     eff  fc-mltime     eff
               run_1:         0.01  100.00%         0.06  100.00%         0.03  100.00%
              run_64:         0.01   54.78%         0.07   96.34%         0.04   75.02%
             run_216:         0.01   42.66%         0.07   96.97%         0.05   66.00%
             run_512:         0.02   37.98%         0.07   93.78%         0.05   58.91%


FAQ:

 * if you see the error
          tt = map(float,r[-1].split())
          ValueError: invalid literal for float(): sh:

   then put removeLines.sh in a directory that is on your path.  Since arcadia cd's around, PWD is not enough.

 * if you see the message "running now?"

   make sure each directory run_1, run_8, etc. has a file screen.out.12345.  (It can be empty.)
