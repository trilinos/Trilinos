% MetaUML, a MetaPost library for typesetting exquisite UML diagrams.
% Copyright (C) 2005 Radu-George Radulescu
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

if known _metauml_instance_mp:
  expandafter endinput
fi;
_metauml_instance_mp:=1;

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

inputonce util_picture;
inputonce util_picture_stack;
inputonce util_shade;

vardef InstanceInfo@#=
  attributes(@#);
  var(color)   foreColor, borderColor;

  FontInfo.@#nameFont(metauml_defaultFont, defaultscale);

  FontInfo.@#attributeFont(metauml_defaultFont, defaultscale);

  ShadeInfo.@#iShade;

  @#foreColor := .9white;
  @#borderColor := black;

  PictureInfo.@#iName      (2, 2, 1, 3)(@#nameFont);

  PictureInfo.@#iAttribute (2, 2, 1.25, 0)(@#attributeFont);

  PictureStackInfo.@#iAttributeStack (2, 2, 2.5, 2)(10.5)(@#iAttribute); % 8.5
  @#iAttributeStack.iPict.bottom := 2;

  @#iName.ignoreNegativeBase := 1;
  @#iAttributeStack.iPict.ignoreNegativeBase := 1;
enddef;

vardef InstanceInfoCopy@#(text src)=
  log "InstanceInfoCopy: Copying Instance";

  attributes(@#);
  var(color)   foreColor, borderColor;

  FontInfoCopy.@#nameFont(src.nameFont);

  FontInfoCopy.@#attributeFont(src.attributeFont);

  ShadeInfoCopy.@#iShade(src.iShade);

  @#foreColor := src.foreColor;
  @#borderColor := src.borderColor;

  PictureInfoCopy.@#iName      (src.iName);

  PictureInfoCopy.@#iAttribute (src.iAttribute);
  PictureStackInfoCopy.@#iAttributeStack (src.iAttributeStack);

enddef;

InstanceInfo.iInstance;

%
% Instance
%
vardef defInstance@#(expr pname) =
  ObjectEquations(@#);
  @#className := "Instance";

  string @#name;
  @#name = pname;

  string @#attributes[];
  numeric @#nAttrs;
  @#nAttrs   := 0;
enddef;

vardef Instance_addAttribute@#(expr pcontent) =

  @#.attributes[@#.nAttrs] := pcontent;
  @#.nAttrs := @#.nAttrs + 1;

enddef;

vardef EInstance@#(text _info)(expr name)(text attributes)=
  log "EInstance begin";
  defInstance@#(name);

  log "Copying Instance info";
  InstanceInfoCopy.@#info(_info);

  for a = attributes:
    log "Adding attribute ", a;
    Instance_addAttribute@#(a);
  endfor;
enddef;

vardef Instance@#(expr name)(text attributes)=
  EInstance@#(iInstance)(name)(attributes);
enddef;

vardef Instance_border@#=
  objectBox(@#)
enddef;

vardef Instance_layout@# =
  if @#layedout = 1:
    log "Instance " & (str @#) & " has already been layed out";
  else:
    @#layedout := 1;
    log "Instance layout: " & (str @#);

    @#maxFeatureWidth := 0;

    EPicture.@#namePict(@#info.iName)(@#name);
    EPictureStack.@#attributeStack(@#info.iAttributeStack)
       (scantokens listArray(@#attributes)(@#nAttrs))("vleftbase");

    layoutObjects(@#namePict, @#attributeStack);

    @#maxFeatureWidth := lmax(@#namePict.width, @#attributeStack.width) + 1;

    @#namePict.midx = @#midx;

    @#namePict.top = @#top - 1;
    @#attributeStack.top = @#namePict.bottom - 1;

    @#attributeStack.left = @#left;

    @#width = @#maxFeatureWidth;
    @#bottom = @#attributeStack.bottom;

    log "Layout done...";
  fi;
enddef;

vardef Instance_drawFeatures@#=
  log "Drawing name";
  drawObject(@#namePict);
  log "Underlining name";
  pair A, B;
  A := @#namePict.se + (-@#namePict.info.left, @#namePict.info.bottom-1);
  B := @#namePict.sw + (@#namePict.info.right, @#namePict.info.bottom-1);
  draw A--B;
  log "Drawing attribute stack";
  drawObject(@#attributeStack);

enddef;

vardef Instance_paintSkin@# =
  log "Painting Instance skin...";

 drawObjectShade(@#);

 fill Instance_border.@# withcolor @#info.foreColor;
 draw Instance_border.@# withcolor @#info.borderColor;

 draw (xpart @#nw, @#attributeStack.top)--(xpart @#se, @#attributeStack.top) withcolor @#info.borderColor;
enddef;

vardef Instance_draw@#=
  log "draw Instance begin...";
  Instance_layout.@#;
  objectEnsurePositioning.@#;
  Instance_paintSkin.@#;
  Instance_drawFeatures.@#;
enddef;

vardef Instance_setDebugMode@#=
  @#.info.iName.boxed := 1;
  @#.info.iAttributeStack.boxed := 1;
  @#.info.iAttributeStack.iPict.boxed := 1;
enddef;
