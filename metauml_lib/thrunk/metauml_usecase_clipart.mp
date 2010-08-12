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

%%
%% This file contains various drawings which can be used in
%% UML usecase diagrams
%%

if known _metauml_usecase_clipart_mp:
  expandafter endinput
fi;
_metauml_usecase_clipart_mp:=1;

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

vardef HumanInfo@#(expr _width, _height, _yHeadRatio, _yBodyRatio, _yHandsRatio, _xHandsRatio, _xLegsRatio)=
  attributes(@#);

  var(numeric) width, height, yHeadRatio, yBodyRatio,
                           yHandsRatio, xHandsRatio, xLegsRatip;
  var(color) foreColor, shadeColor;
  var(numeric) shade;

  @#width := _width;
  @#height := _height;
  @#yHeadRatio := _yHeadRatio;
  @#yBodyRatio := _yBodyRatio;
  @#yHandsRatio := _yHandsRatio;
  @#xHandsRatio := _xHandsRatio;
  @#xLegsRatio := _xLegsRatio;

  @#foreColor := black;
  @#shadeColor := .7white;
  @#shade := .5;

  Margins.@#(2,2,1,0);
enddef;

vardef HumanInfoCopy@#(text src)=
  HumanInfo@#(src.width, src.height, src.yHeadRatio, src.yBodyRatio, src.yHandsRatio, src.xHandsRatio, src.xLegsRatio);

  _n_ := str @#;
  var(color) foreColor, shadeColor;
  var(numeric) shade;

  @#foreColor := src.foreColor;
  @#shadeColor := src.shadeColor;
  @#shade := src.shade;

  MarginsCopy.@#(src);
enddef;

HumanInfo.iHuman(25, 25, .3, .35, .55, 1, 1);

vardef EHuman@#(text _info)=
  ObjectEquations(@#);
  @#className := "Human";

  var(path) pHead, pBody, pHands, pLegs;

  HumanInfoCopy.@#info(_info);
enddef;

vardef Human@#=
  EHuman@#(iHuman);
enddef;

vardef Human_layout@#=
  if @#layedout = 1:
    log "Human " & (str @#) & " has already been layed out";
  else:
    @#layedout := 1;
    log "Human layout: " & (str @#);

    attributes(@#);

    @#height = @#info.height;
    @#width  = @#info.width;

    log "(W, H) MUST BE SET. Their values are: " & (decimal @#height) &" "& (decimal @#width);

    var(numeric) actualHeight, actualWidth;
    @#actualHeight = @#height - @#info.top - @#info.bottom;
    @#actualWidth = @#width - @#info.left - @#info.right;

    % HEAD
    var(pair) headS, headE, headW;

    ypart @#headW = ypart @#headE;
    .5[@#n - (0, @#info.top), @#headS] = .5[@#headE, @#headW];
    (ypart @#n - @#info.top) - (ypart @#headS) =
        (xpart @#headE) - (xpart @#headW) =
        @#actualHeight * @#info.yHeadRatio;
    (xpart @#headS) = (xpart @#n);

    % BODY
    var(pair) bodyS;
    xpart @#bodyS = xpart @#headS;
    ypart @#headS - ypart @#bodyS = @#actualHeight * @#info.yBodyRatio;

    % HANDS
    var(pair) handsW, handsE;
    ypart @#handsW = ypart @#handsE = ypart @#n - (@#actualHeight * @#info.yHandsRatio);
    .5[xpart @#handsE, xpart @#handsW] = @#midx;
    xpart @#handsE - xpart @#handsW = @#actualWidth * @#info.xHandsRatio;

    % LEGS
    var(pair) footW, footE;
    ypart @#footW = ypart @#footE = ypart @#s + @#info.bottom;
    .5[xpart @#footW, xpart @#footE] = @#midx;
    xpart @#footE - xpart @#footW = @#actualWidth * @#info.xLegsRatio;

    log "Human layout completed";
  fi;
enddef;

vardef Human_draw@#=
  log "Drawing human";

  Human_layout.@#;
  objectEnsurePositioning.@#;

  % all coordinates must be set to known values in order to compute the paths
  @#pHead := (@#n - (0, @#info.top)) .. @#headE .. @#headS .. @#headW .. cycle;
  @#pBody := @#headS -- @#bodyS;
  @#pHands := @#handsW -- @#handsE;
  @#pLegs := @#footW -- @#bodyS -- @#footE;

  attributes(@#);
  var(pair) shadeShift;
  @#shadeShift = (@#info.shade, -@#info.shade);

  draw @#pHead shifted @#shadeShift withcolor @#info.shadeColor;
  draw @#pBody shifted @#shadeShift withcolor @#info.shadeColor;
  draw @#pHands shifted @#shadeShift withcolor @#info.shadeColor;
  draw @#pLegs shifted @#shadeShift withcolor @#info.shadeColor;

  draw @#pHead withcolor @#info.foreColor;
  draw @#pBody withcolor @#info.foreColor;
  draw @#pHands withcolor @#info.foreColor;
  draw @#pLegs withcolor @#info.foreColor;
enddef;

vardef Human_border@#=
  objectBox(@#)
enddef;
