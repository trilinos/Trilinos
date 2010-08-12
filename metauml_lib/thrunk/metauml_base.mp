% MetaUML, a MetaPost library for typesetting exquisite UML diagrams.
% Copyright (C) 2006 Ovidiu Gheorghies, Radu-George Radulescu
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

if known _metauml_base_mp:
  expandafter endinput
fi;
_metauml_base_mp:=1;

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

vardef markPoint(text p)(text color) =
	draw p withpen pencircle scaled 1 withcolor color;
enddef;

vardef drawLine(text pA)(text pB)(text base) =
	draw pA--pB;
enddef;

vardef drawArrow(text pA)(text pB)(text base) =
  numeric dx, dy;
  dx = xpart(pB - pA);
  dy = ypart(pB - pA);

  numeric alfa;
  alfa = angle(dx, dy);

  pair pC, pD;
  pC := (0, base/2);
  pC := pC rotated (alfa) shifted pA;
  pD := (0, base/2);
  pD := pD rotated (180+alfa) shifted pA;

  %fill pA--pC--pB--pD--cycle withcolor .8white;

  draw pA--pB;
  draw pB--pC;
  draw pB--pD;
enddef;

vardef drawTriangle(text pA)(text pB)(text base) =
  numeric dx, dy;
  dx = xpart(pB - pA);
  dy = ypart(pB - pA);

  numeric alfa;
  alfa = angle(dx, dy);
  pair pC, pD;
  pC := (0, base/2);
  pC := pC rotated (alfa) shifted pA;
  pD := (0, base/2);
  pD := pD rotated (180+alfa) shifted pA;

  fill pA--pC--pB--pD--cycle withcolor white;

  draw pA--pC;% withcolor red;
  draw pA--pD;% withcolor green;
  draw pB--pC;% withcolor blue;
  draw pB--pD;
enddef;

vardef drawDiamond(text pA)(text pB)(text height) =
  numeric dx, dy;
  dx = xpart(pB - pA);
  dy = ypart(pB - pA);

  numeric alfa;
  alfa = angle(dx, dy);

  numeric width;
  width = _length(pA)(pB);

  path diamond;
  diamond = (0,0)--(width/2, height/2)--(width, 0)--(width/2, -height/2)--cycle;

  fill diamond rotated alfa shifted pA withcolor white;
  draw diamond rotated alfa shifted pA;
enddef;

vardef drawDiamondBlack(text pA)(text pB)(text height) =
  numeric dx, dy;
  dx = xpart(pB - pA);
  dy = ypart(pB - pA);

  numeric alfa;
  alfa = angle(dx, dy);

  numeric width;
  width = _length(pA)(pB);

  path diamond;
  diamond = (0,0)--(width/2, height/2)--(width, 0)--(width/2, -height/2)--cycle;

  fill diamond rotated alfa shifted pA;
enddef;

vardef drawCrossedCircle(text pA)(text pB)(text height) =
  numeric dx, dy;
  dx = xpart(pB - pA);
  dy = ypart(pB - pA);

  numeric alfa;
  alfa = angle(dx, dy);

  numeric width;
  width = _length(pA)(pB);

  pair pC, pD;
  pC := (width/2, height/2 + 1);
  pC := pC rotated (alfa) shifted pA;
  pD := (width/2, -height/2 - 1);
  pD := pD rotated (alfa) shifted pA;

  fill pA..pB..pA -- cycle withcolor white;

  draw pA -- pB;
  draw pC -- pD;
  draw pA .. pB .. pA;
enddef;

vardef drawSemiCircle(text pA)(text pB)(text height) =
  numeric dx, dy;
  dx = xpart(pB - pA);
  dy = ypart(pB - pA);

  numeric alfa;
  alfa = angle(dx, dy);

  numeric width;
  width = _length(pA)(pB);

  pair pC, pD;
  pC := (width/2 + 5, height/2 + 1);
  pC := pC rotated (alfa) shifted pB;
  pD := (width/2 + 5, -height/2 - 1);
  pD := pD rotated (alfa) shifted pB;

  draw pC .. pA .. pD;
enddef;

vardef drawCircle(text pA)(text pB)(text height) =
  fill pA..pB..pA -- cycle withcolor white;

  draw pA .. pB .. pA;
enddef;

vardef drawLine(expr my_path) =
  draw my_path;
enddef;

vardef drawLineDashed(expr my_path) =
  draw my_path dashed evenly;
enddef;

vardef drawNothing(text pA)(text pB)(text base) =
enddef;

def shorterLineEndPoint(text pA)(text pB)(text delta) =
  (xpart(pA) +  _dx(pA)(pB)*(1 - delta/_length(pA)(pB)), ypart(pA) + _dy(pA)(pB) * (1 - delta/_length(pA)(pB)) )
enddef;
