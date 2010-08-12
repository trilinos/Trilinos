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

if known _metauml_behavioral_common_mp:
  expandafter endinput
fi;
_metauml_behavioral_common_mp:=1;

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

vardef SpecialStateInfo@#=
  numeric @#diameter;
  color @#discColor;

  numeric @#shade;
  color @#shadeColor, @#foreColor;

  @#diameter := 0.4cm;
  @#discColor := black;
  @#foreColor := black;

  ShadeInfo.@#iShade;
enddef;

vardef SpecialStateInfoCopy@#(text src)=
  SpecialStateInfo.@#;
  ShadeInfo_copy.@#iShade(src.iShade);

  @#diameter := src.diameter;
  @#discColor := src.discColor;
  @#foreColor := src.foreColor;
enddef;

%%%%%%%%%%%%%%%%%%
%     BEGIN      %
%%%%%%%%%%%%%%%%%%

SpecialStateInfo.iBegin;

vardef Begin@#=
  ObjectEquations(@#);
  @#className := "Begin";

  SpecialStateInfoCopy.@#info(iBegin);
enddef;

vardef Begin_layout@#=
  if @#layedout = 1:
    log "Begin " & (str @#) & " has already been layed out";
  else:
    @#layedout := 1;

    @#width = @#info.diameter;
    @#height = @#info.diameter;
  fi;
enddef;

vardef Begin_draw@#=
  Begin_layout@#;
  objectEnsurePositioning.@#;

  path @#border;
  @#border := @#n..@#e..@#s..@#w..cycle;

  drawObjectShade(@#);

  fill @#border withcolor @#info.discColor;
  draw @#border withcolor @#info.foreColor;
enddef;

vardef Begin_border@#=
  @#border
enddef;

%%%%%%%%%%%%%%%%%%
%      END       %
%%%%%%%%%%%%%%%%%%

SpecialStateInfo.iEnd;
iEnd.innerSpacing := 2;
iEnd.discColor := white;

vardef End@#=
  ObjectEquations(@#);
  @#className := "End";
  SpecialStateInfoCopy.@#info(iEnd);
  @#info.innerSpacing := iEnd.innerSpacing;
enddef;

vardef End_layout@#=
  Begin_layout@#;
enddef;

vardef End_draw@#=
  Begin_draw.@#;

  path @#innerDisc;
  @#innerDisc := @#n - (0, @#info.innerSpacing) ..
                 @#e - (@#info.innerSpacing, 0) ..
                 @#s + (0, @#info.innerSpacing) ..
		 @#w + (@#info.innerSpacing, 0) ..cycle;

  fill @#innerDisc withcolor @#info.foreColor;
enddef;

vardef End_border@#=
  @#border
enddef;

%%%%%%%%%%%%%%%%%%
%   ENTRY POINT  %
%%%%%%%%%%%%%%%%%%

SpecialStateInfo.iEntryPoint;
iEntryPoint.discColor := white;
iEntryPoint.innerSpacing := 2;

vardef EntryPoint@#=
  ObjectEquations(@#);
  @#className := "EntryPoint";
  SpecialStateInfoCopy.@#info(iEntryPoint);
  @#info.innerSpacing := iEntryPoint.innerSpacing;
enddef;

vardef EntryPoint_layout@#=
  Begin_layout@#;
enddef;

vardef EntryPoint_draw@#=
  Begin_draw.@#;
enddef;

vardef EntryPoint_border@#=
  @#border
enddef;

%%%%%%%%%%%%%%%%%%
%   EXIT POINT   %
%%%%%%%%%%%%%%%%%%

SpecialStateInfo.iExitPoint;
iExitPoint.discColor := white;
iExitPoint.innerSpacing := 2;

vardef ExitPoint@#=
  ObjectEquations(@#);
  @#className := "ExitPoint";
  SpecialStateInfoCopy.@#info(iExitPoint);
  @#info.innerSpacing := iExitPoint.innerSpacing;
  @#info.iShade.shift := 0;
enddef;

vardef ExitPoint_layout@#=
  Begin_layout@#;
enddef;

vardef ExitPoint_draw@#=
  Begin_draw.@#;

  _n_ := str @#;
  generic_declare(pair) _n.urt, _n.lrt, _n.llft, _n.ulft;

  @#urt  := @#c + @#info.diameter/2 * (sind(45), cosd(45));
  @#lrt  := @#c + @#info.diameter/2 * (sind(135), cosd(135));
  @#llft := @#c + @#info.diameter/2 * (sind(225), cosd(225));
  @#ulft := @#c + @#info.diameter/2 * (sind(315), cosd(315));

  draw @#urt -- @#llft;
  draw @#ulft -- @#lrt;
enddef;

vardef ExitPoint_border@#=
  @#border
enddef;

%%%%%%%%%%%%%%%%%%
%   FlowFinal    %
%%%%%%%%%%%%%%%%%%

SpecialStateInfo.iFlowFinal;
iFlowFinal.discColor := white;
iFlowFinal.innerSpacing := 2;

vardef FlowFinal@#=
  ObjectEquations(@#);
  @#className := "FlowFinal";
  SpecialStateInfoCopy.@#info(iFlowFinal);
  @#info.innerSpacing := iFlowFinal.innerSpacing;
  %@#info.iShade.shift := 0;
enddef;

vardef FlowFinal_layout@#=
  Begin_layout@#;
enddef;

vardef FlowFinal_draw@#=
    ExitPoint_draw@#
enddef;

vardef FlowFinal_border@#=
  @#border
enddef;
