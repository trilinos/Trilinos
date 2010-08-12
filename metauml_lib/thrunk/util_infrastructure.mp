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

if known _util_infrastructure_mp:
  expandafter endinput
fi;
_util_infrastructure_mp:=1;

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

vardef attributes(text objectName)=
  _n_ := str objectName
enddef;

vardef variables_ext text variables=
  save ret, firstVar;
  string ret;
  ret := "";

  numeric firstVar;
  firstVar := 1;

  forsuffixes v = variables:
    if firstVar = 0:
      ret := ret & ",";
    else:
      firstVar := 0;
    fi;

    ret := ret & " _n." & (str v);
  endfor;

  ret
enddef;

vardef var(text type) text variables =
  save _variables;
  string _variables;
  _variables := variables_ext variables;

  generic_declare(type) scantokens _variables;
enddef;


