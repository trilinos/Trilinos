% MetaUML, a MetaPost library for typesetting exquisite UML diagrams.
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

if known _metauml_class_clipart_mp:
  expandafter endinput
fi;
_metauml_class_clipart_mp:=1;

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

%%
%% This file contains various drawing which can be used in
%% UML class diagrams
%%

vardef LockInfo@#(expr _width, _height, _y, _d, _span, _faceColor, _topColor, _rightColor, _lockColor)=
  numeric @#y, @#d, @#span, @#width, @#height;
  color @#faceColor, @#topColor, @#rightColor, @#lockColor;

  @#width := _width;
  @#height := _height;

  @#heightRatio := _y;
  @#depthRatio := _d;
  @#span := _span;

  @#faceColor := _faceColor;
  @#topColor := _topColor;
  @#rightColor := _rightColor;
  @#lockColor := _lockColor;
enddef;

vardef LockInfoCopy@#(text src)=
  LockInfo@#(src.width, src.height, src.heightRatio, src.depthRatio, src.span,
    src.faceColor, src.topColor, src.rightColor, src.lockColor);
enddef;

LockInfo.iLock(8, 8, .55, .15, .55, .4white, .6white, .7white, .3white);

def lockLockPen=pencircle scaled .8pt enddef;
def lockLockTension=tension(1) enddef;

vardef ELock@#(text lockInfo)(expr how)=
  ObjectEquations(@#);
  @#className := "Lock";

  generic_declare(path) _n.pFace, _n.pRight, _n.pTop;
  generic_declare(path) _n.pLock;
  generic_declare(string) _n.mode;

  LockInfoCopy.@#info(lockInfo);

  @#height = @#info.height;
  @#width  = @#info.width;

  @#mode := how;

  log "ELock completed";
enddef;

vardef Lock@#(expr how)=
  ELock@#(iLock)(how);
enddef;

vardef Lock_draw@#=
  log "Drawing lock";
  log @#nw;
  log @#se;
  generic_declare(numeric) _n.depth,_n.ya, _n.xa;

  @#depth = @#info.width * @#info.depthRatio;
  @#ya = ypart @#nw - @#height * @#info.heightRatio;
  @#xa = xpart @#se - @#depth;

  @#pFace  = (xpart @#nw, @#ya) -- (@#xa, @#ya) -- (@#xa, ypart @#se) -- (xpart @#nw, ypart @#se) -- cycle;
  @#pRight = (@#xa, ypart @#se) -- (xpart @#se, ypart @#se + @#depth) --
              (xpart @#se, @#ya + @#depth) -- (@#xa, @#ya) -- cycle;
  @#pTop   = (@#xa, @#ya) -- (xpart @#se, @#ya + @#depth) --
              (xpart @#nw + @#depth, @#ya + @#depth) -- (xpart @#nw, @#ya) -- cycle;

  save myMode; string myMode;
  myMode = "base";
  if ((@#mode="private") or (@#mode="-")):
    myMode := "closed";
  elseif ((@#mode="protected") or (@#mode="#")):
    myMode := "opened";
  fi;

  if myMode="closed":
    log "Closed!";
    generic_declare(pair) _n.A, _n.B, _n.C, _n.M;
    generic_declare(numeric) _n.span;
    @#span = @#info.span * @#width;

    @#M = ((xpart @#se + xpart @#nw) / 2, @#ya + @#depth/2);
    log @#M;
    @#A = (xpart @#M - (@#span / 2), ypart @#M);
    log @#A;
    @#B = (xpart @#M + (@#span / 2), ypart @#M);
    log @#B;
    @#C = (xpart @#M, ypart @#nw);
    log @#C;
    linecap:=butt;

    @#pLock = @#A {curl 0}..lockLockTension .. @#C .. lockLockTension .. {curl 0} @#B;
         %withcolor @#info.lockColor withpen lockLockPen;
  elseif myMode="opened":
    log "Opened!";

    generic_declare(pair) _n.A, _n.B, _n.C, _n.M;
    generic_declare(numeric) _n.span;
    @#span = @#info.span * @#width;

    @#M = ((xpart @#se + xpart @#nw) / 2, @#ya + @#depth/2);
    @#A = (xpart @#M - (@#span / 2), ypart @#M) + (0, @#height/3);
    @#B = (xpart @#M + (@#span / 2), ypart @#M);
    @#C = (xpart @#M, ypart @#nw);
    log @#M, @#A, @#B, @#C;

    linecap:=butt;

    @#pLock = @#A {curl 0}..lockLockTension .. @#C .. lockLockTension .. {curl 0} @#B;
         %withcolor @#info.lockColor withpen lockLockPen;
  fi;

  fill @#pFace withcolor @#info.faceColor;
  fill @#pRight withcolor @#info.rightColor;
  fill @#pTop withcolor @#info.topColor;

  if not ((myMode = "opened") or (myMode="closed")):
    log "No lock to draw...";
  else:
    draw @#pLock withcolor @#info.lockColor withpen lockLockPen;
  fi;
enddef;


