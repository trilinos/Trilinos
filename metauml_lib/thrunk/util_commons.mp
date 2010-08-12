% MetaUtil, an easier MetaPost life
% Copyright (C) 2005 Ovidiu Gheorghies
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

if known _util_commons_mp:
  expandafter endinput
fi;
_util_commons_mp:=1;

% Sadly, this copy of the macro is needed to prevent multiple file loads being shown by MetaPost.
% The guard values (such as _metauml_mp) do ensure that the file isn't loaded multiple times, 
% but this macro makes sure that MetaPost won't try to load the file and display a message for that.
def inputonce text libraryFile=
	if not known scantokens ("_" & str libraryFile & "_mp"):
		%includeonce% show "Loading " & str libraryFile;
		scantokens ("input " & str libraryFile);
	else:
		%includeonce% show str libraryFile & " already loaded.";
	fi;
enddef;

inputonce util_log;

vardef lmax(text items)=
  log "finding lmax of items " & str items;
  save current, nItems;
  numeric current, nItems;
  nItems := 0;
  current:= 0;

  for item = items:
    log nItems;
    log item;

    if nItems = 0:
      current := item;
    else:
      current := max(current)(item);
    fi;
    nItems := nItems + 1;
  endfor;

  current
enddef;

vardef lmin(text items)=
  save current, nItems;
  numeric current, nItems;
  nItems := 0;
  current:= 0;

  for item = items:
    if nItems = 0:
      current := item;
    else:
      current := min(current)(item);
    fi;
    nItems := nItems + 1;
  endfor;

  current
enddef;

def max(text a)(text b)=
  if (a > b): a%
  else:b%
  fi;
enddef;

def min(text a)(text b)=
  if (a < b): a%
  else:b%
  fi;
enddef;

vardef ensurePict(text pictOrText)(expr font, scale) =
  save p; picture p;
  if picture pictOrText:  p=pictOrText
  else: p = pictOrText infont font scaled scale
  fi;
  p
enddef;

def pictHeight(text obj) =
  (ypart (urcorner obj) - ypart (llcorner obj))
enddef;

def pictWidth(text obj) =
    (xpart (urcorner obj) - xpart (llcorner obj))
enddef;

vardef listArray(text array)(text nElements)=
  save objEnum;
  string objEnum;
  objEnum := "";

  for i = 0 upto nElements-1:
    if i>0:
      objEnum := objEnum & ", ";
    fi;
    objEnum := objEnum & (str array) & (decimal i);
  endfor;

  objEnum
enddef;

vardef enumToString(text enumeration)(expr prefix)=
  save ret, firstVar;
  string ret;
  ret := "";

  numeric firstVar;
  firstVar := 1;

  forsuffixes v = enumeration:
    if firstVar = 0:
      ret := ret & ",";
    else:
      firstVar := 0;
    fi;

    ret := ret & prefix & (str v);
  endfor;

  ret
enddef;

def _dx(text pA)(text pB) =
  (xpart(pB)-xpart(pA))
enddef;

def _dy(text pA, pB) =
  (ypart(pB)-ypart(pA))
enddef;

def _length(text pA)(text pB) =
  sqrt (_dx(pA)(pB)*_dx(pA)(pB) + _dy(pA)(pB)*_dy(pA)(pB))
enddef;
