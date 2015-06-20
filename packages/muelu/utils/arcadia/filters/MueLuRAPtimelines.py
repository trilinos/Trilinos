#!/usr/bin/env python
#MueLu MM coarse granulariy timelines
LABELS       = ['RAP total', 'AxP', 'Rx(AP)']                      # analysis header labels
TIMELINES    = ['MueLu: RAPFactory : Computing Ac', 'MueLu: RAPFactory: MxM: A x P (sub, total)', 'MueLu: RAPFactory: MxM: R x (AP) (explicit) (sub, total)']
def PARSEFUNC(fileName,stringToFind):
  # sed kungfu to get rid of bad parentheses, since cut uses parens to find the timers
  return "grep -i \"" + stringToFind + "\" " + fileName \
             + "| tail -n 1" \
             + "| sed \"s/(total)/[total]/\"" \
             + "| sed \"s/(explicit)/[explicit]/\"" \
             + "| sed \"s/(sub, total)/[sub, total]/\"" \
             + "| sed \"s/(sub/[sub/\"" \
             + "| sed \"s/(AP)/AP/\"" \
             + "| sed \"s/\(level=.\))/\\1]/\"" \
             + "| cut -f3 -d')' | cut -f1 -d'('"
