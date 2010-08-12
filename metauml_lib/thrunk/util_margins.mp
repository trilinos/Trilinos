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

if known _util_margins_mp:
   expandafter endinput
fi;
_util_margins_mp:=1;

vardef Margins@#(expr pleft, pright, ptop, pbottom)=
  attributes(@#);
  var(numeric) left, right, top, bottom;

  @#left := pleft;
  @#right := pright;
  @#top := ptop;
  @#bottom := pbottom;
enddef;

vardef MarginsCopy@#(text src)=
  Margins@#(src.left, src.right, src.top, src.bottom);
enddef;
