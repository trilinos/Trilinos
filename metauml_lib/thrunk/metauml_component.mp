% MetaUML, a MetaPost library for typesetting exquisite UML diagrams.
% Copyright (C) 2006 Radu-George Radulescu, Ovidiu Gheorghies
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

if known _metauml_component_mp:
  expandafter endinput
fi;
_metauml_component_mp:=1;

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

inputonce metauml_defaults;
inputonce util_log;

vardef ComponentInfo@#=
  attributes(@#);
  var(color)   foreColor, borderColor;
  var(numeric) minimumNameContentsDifference;

  FontInfo.@#iFont(metauml_defaultFont, defaultscale);
  @#iFont.ignoreNegativeBase := 1;

  PictureInfo.@#iText  (2, 2, 2, 2)(@#iFont);
  PictureStackInfo.@#iNameStack (2, 2, 2, 2)(9)(@#iText);
  @#iNameStack.iPict.ignoreNegativeBase := 1;

  Margins.@#(2,2,2,2);

  ShadeInfo.@#iShade;
  @#foreColor := .9white;
  @#borderColor := black;

  @#minimumNameContentsDifference := 20;

  GroupInfo.@#iContentsGroup(2, 2, 10, 2);
enddef;

vardef ComponentInfoCopy@#(text src)=
  log "ComponentInfoCopy: Copying component info";

  attributes(@#);
  var(color)   foreColor, borderColor;
  var(numeric) minimumNameContentsDifference;

  PictureStackInfoCopy.@#iNameStack (src.iNameStack);
  MarginsCopy.@#(src);

  ShadeInfoCopy.@#iShade(src.iShade);
  @#foreColor := src.foreColor;
  @#borderColor := src.borderColor;

  @#minimumNameContentsDifference := src.minimumNameContentsDifference;

  GroupInfoCopy.@#iContentsGroup(src.iContentsGroup);
enddef;

ComponentInfo.iComponent;

vardef Component@#(text contents)(text _subItems)=
  EComponent.@#(iComponent)(contents)(_subItems);
enddef;

vardef EComponent@#(text _info)(text contents)(text _subItems)=
  ObjectEquations(@#);
  @#className := "Component";

  ComponentInfoCopy.@#info(_info);

  numeric @#nLines; @#nLines := 0;
  string @#lines[];
  numeric @#nSubItems; @#nSubItems := 0;

  for l=contents:
    @#lines[@#nLines] := l;
    @#nLines := @#nLines + 1;
  endfor;

  for s=_subItems:
    @#nSubItems := @#nSubItems + 1;
  endfor;

  EGroup.@#subItems(@#info.iContentsGroup)(_subItems);

  EPictureStack.@#nameStack(@#info.iNameStack)(scantokens listArray(@#lines)(@#nLines))("vcenterbase");
enddef;

vardef Component_layout@#=
  if @#layedout = 1:
    log "Component '" & (str @#) & "' has already been layed out";
  else:
    @#layedout := 1;
    log "Component layout: '" & (str @#) & "'";

    log "Component name layout '" & (str @#) & "'...";
    PictureStack_layout.@#nameStack;
    log "SubItems for component layout '" & (str @#) & "'...";
    Group_layout.@#subItems;

    log "All elements in component '" & (str @#) & "' successfully layed out, integrating...";


    if @#subItems.width < @#nameStack.width + @#info.minimumNameContentsDifference:
      @#contentWidth = @#nameStack.width + @#info.minimumNameContentsDifference;
    else:
      @#contentWidth = lmax(@#nameStack.width, @#subItems.width);
    fi;
    
    @#width = @#contentWidth + @#info.left + @#info.right + 15;
    
    if @#nSubItems = 0:
       @#height = @#info.top + @#info.bottom + @#subItems.height + 15;
       @#nameStack.n = @#n - (0, 7);
    else:
       @#height = @#info.top + @#info.bottom + @#subItems.height + 30;
       @#nameStack.n = @#n + (0, @#info.top) - (0, 5);
    fi;

    @#subItems.n = @#n + (0, -((@#height - @#subItems.height)/2)) - (0, 5);

    log "Component layout for " & (str @#) & " ended.";
  fi;
enddef;

vardef Component_draw@#=
  Component_layout.@#;
  objectEnsurePositioning.@#;

  attributes(@#);
  var(path) border;

  @#border := @#ne -- @#nw -- @#sw -- @#se -- cycle;

  drawObjectShade(@#);

  fill @#border withcolor @#info.foreColor;
  draw @#border withcolor @#info.borderColor;
  
  drawComponentVisualStereotype(@#ne);
  
  drawObject(@#nameStack);
  drawObject(@#subItems);

enddef;

vardef drawComponentVisualStereotype(text ne)=
  % lud = left-up-down of the top-right rectangle; luu = left-up-up of the top-right rectangle etc.
  pair lud, luu, ru, rd, ldd, ldu, lmd, lmu;
  % urectld = upper rectangle left-down, lrectru = lower rectangle right-up
  pair urectld, urectlu, urectru, urectrd, lrectld, lrectlu, lrectru, lrectrd;
  lud = ne - (14, 5 );
  luu = lud  + (0,  2 );
  ru  = luu  + (8,  0 );
  rd  = ru   - (0,  14);
  ldd = rd   - (8,  0 );
  ldu = ldd  + (0,  2 );
  lmd = ldu  + (0,  3 );
  lmu = lud  - (0,  3 );

  urectld = lmu - (3, 0);
  urectlu = lud - (3, 0);
  urectru = lud + (3, 0);
  urectrd = lmu + (3, 0);

  lrectld = ldu - (3, 0);
  lrectlu = lmd - (3, 0);
  lrectru = lmd + (3, 0);
  lrectrd = ldu + (3, 0);

  draw lud -- luu -- ru -- rd -- ldd -- ldu;
  draw lmd -- lmu;
  draw urectld -- urectlu -- urectru -- urectrd -- urectld;
  draw lrectld -- lrectlu -- lrectru -- lrectrd -- lrectld;
enddef;

vardef Component_border@#=
  @#border
enddef;
