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

if known _metauml_package_mp:
  expandafter endinput
fi;
_metauml_package_mp:=1;

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

vardef PackageInfo@#=
  attributes(@#);
  var(color)   foreColor, borderColor;
  var(numeric) forceEmptyContent;
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

  @#forceEmptyContent := 0;

  @#minimumNameContentsDifference := 20;

  GroupInfo.@#iContentsGroup(2, 2, 10, 2);
enddef;

vardef PackageInfoCopy@#(text src)=
  log "PackageInfoCopy: Copying package info";

  attributes(@#);
  var(color)   foreColor, borderColor;
  var(numeric) forceEmptyContent;
  var(numeric) minimumNameContentsDifference;

  PictureStackInfoCopy.@#iNameStack (src.iNameStack);
  MarginsCopy.@#(src);

  ShadeInfoCopy.@#iShade(src.iShade);
  @#foreColor := src.foreColor;
  @#borderColor := src.borderColor;

  @#forceEmptyContent := src.forceEmptyContent;

  @#minimumNameContentsDifference := src.minimumNameContentsDifference;

  GroupInfoCopy.@#iContentsGroup(src.iContentsGroup);
enddef;

PackageInfo.iPackage;

vardef Package@#(text contents)(text _subItems)=
  EPackage.@#(iPackage)(contents)(_subItems);
enddef;

vardef EPackage@#(text _info)(text contents)(text _subItems)=
  ObjectEquations(@#);
  @#className := "Package";

  PackageInfoCopy.@#info(_info);

  numeric @#nLines; @#nLines := 0;
  string @#lines[];
%  numeric @#nSubItems; @#nSubItems := 0;
  numeric @#nameInCenter; @#nameInCenter := 0;

  for l=contents:
    @#lines[@#nLines] := l;
    @#nLines := @#nLines + 1;
  endfor;

%  for s=_subItems:
%    @#nSubItems := @#nSubItems + 1;
%  endfor;

  EGroup.@#subItems(@#info.iContentsGroup)(_subItems);

  EPictureStack.@#nameStack(@#info.iNameStack)(scantokens listArray(@#lines)(@#nLines))("vcenterbase");
enddef;

vardef Package_layout@#=
  if @#layedout = 1:
    log "Package '" & (str @#) & "' has already been layed out";
  else:
    @#layedout := 1;
    log "Package layout: '" & (str @#) & "'";

    log "Package name layout '" & (str @#) & "'...";
    PictureStack_layout.@#nameStack;
    log "SubItems for package layout '" & (str @#) & "'...";
    Group_layout.@#subItems;

    log "All elements in package '" & (str @#) & "' successfully layed out, integrating...";

    if @#subItems.nObjects = 0:
        if @#info.forceEmptyContent = 0:
          @#nameInCenter := 1;
        fi;
    fi;

    if @#subItems.width < @#nameStack.width + @#info.minimumNameContentsDifference:
      @#contentWidth = @#nameStack.width + @#info.minimumNameContentsDifference;
    else:
      @#contentWidth = lmax(@#nameStack.width, @#subItems.width);
    fi;

    if @#nameInCenter = 1:
      @#height = @#info.top + @#info.bottom + @#nameStack.height;
      @#width = @#info.left + @#info.right + @#nameStack.width;
      @#nameStack.sw = @#sw + (@#info.left, @#info.bottom);
    else:
      @#nameStack.sw = @#nw;
      @#width = @#contentWidth + @#info.left + @#info.right;
      @#height = @#info.top + @#info.bottom + @#subItems.height;
    fi;

    @#subItems.n = @#n + (0, -((@#height - @#subItems.height)/2));

    log "Package layout for " & (str @#) & " ended.";
  fi;
enddef;

vardef Package_draw@#=
  Package_layout.@#;
  objectEnsurePositioning.@#;

  attributes(@#);
  var(path) border, nameStackBorder, emptyNameBorder;

  @#border := @#ne -- @#nw -- @#sw -- @#se -- cycle;
  @#nameStackBorder := @#nameStack.ne -- @#nameStack.nw -- @#nameStack.sw -- @#nameStack.se -- cycle;
  @#emptyNameBorder := @#nw -- @#nw + (0, 10) -- @#n + (0, 10) -- @#n -- cycle;

  drawObjectShade(@#);

  fill @#border withcolor @#info.foreColor;
  draw @#border withcolor @#info.borderColor;
  if @#nameInCenter = 0:
     fill @#nameStackBorder withcolor @#info.foreColor;
     draw @#nameStackBorder withcolor @#info.borderColor;
  else:
     fill @#emptyNameBorder withcolor @#info.foreColor;
     draw @#emptyNameBorder withcolor @#info.borderColor;
  fi;


  drawObject(@#nameStack);
  drawObject(@#subItems);

enddef;

vardef Package_border@#=
  @#border
enddef;
