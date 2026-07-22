
/*
 * Copyright (C) 1991, 1992, 1993, 2022, 2023, 2024, 2024 by Chris Thewalt (thewalt@ce.berkeley.edu)
 *
 * Permission to use, copy, modify, and distribute this software
 * for any purpose and without fee is hereby granted, provided
 * that the above copyright notices appear in all copies and that both the
 * copyright notice and this permission notice appear in supporting
 * documentation.  This software is provided "as is" without express or
 * implied warranty.
 *
 * Thanks to the following people who have provided enhancements and fixes:
 *   Ron Ueberschaer, Christoph Keller, Scott Schwartz, Steven List,
 *   DaviD W. Sanderson, Goran Bostrom, Michael Gleason, Glenn Kasten,
 *   Edin Hodzic, Eric J Bivona, Kai Uwe Rommel, Danny Quah, Ulrich Betzler
 */

#pragma once

namespace SEAMS {
  char *getline_int(const char *); /* read a line of input */
  void  gl_setwidth(size_t);       /* specify width of screen */
  void  gl_histadd(const char *);  /* adds entries to hist */
} // namespace SEAMS
