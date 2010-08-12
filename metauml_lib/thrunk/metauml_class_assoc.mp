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

if known _metauml_class_assoc_mp:
  expandafter endinput
fi;
_metauml_class_assoc_mp:=1;

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
inputonce metauml_defaults;

FontInfo.assocFont(metauml_defaultFont, .7);
PictureInfo.iAssoc(2, 2, 2, 2)(assocFont);
iAssoc.ignoreNegativeBase := 1;

vardef QualifiedAssociationInfo@#=
  color @#background;
  color @#borderColor;
  numeric @#borderScale;

  @#background := .9white;
  @#borderColor := black;
  @#borderScale := .5;

  FontInfo.@#iFont(metauml_defaultFont, .9);
  PictureInfo.@#iPict(2, 2, 2, 2)(@#iFont);
  PictureStackInfo.@#iStack(2, 2, 1, 1)(7)(@#iPict);

  @#iStack.iPict.ignoreNegativeBase := 1;
enddef;

vardef QualifiedAssociationInfoCopy@#(text src)=
  color @#background;
  color @#borderColor;
  numeric @#borderScale;

  @#background := src.background;
  @#borderColor := src.borderColor;
  @#borderScale := src.borderScale;

  FontInfoCopy.@#iFont(src.iFont);
  PictureInfoCopy.@#iPict(src.iPict);
  PictureStackInfoCopy.@#iStack(src.iStack);
enddef;

QualifiedAssociationInfo.iQualifiedAssociation;

%
% QualifiedAssociation
%
vardef EQualifiedAssociation@#(text _info)(text qualifiedAssociations) =
  ObjectEquations(@#);
  @#className := "QualifiedAssociation";

  QualifiedAssociationInfoCopy.@#info(_info);

  EPictureStack.@#elementsStack(@#info.iStack)(qualifiedAssociations)("vleftbase");

  @#nw = @#elementsStack.nw;
  @#se = @#elementsStack.se;
enddef;

vardef QualifiedAssociation@#(text qualifiedAssociations)=
  EQualifiedAssociation@#(iQualifiedAssociation)(qualifiedAssociations);
enddef;

vardef QualifiedAssociation_paintSkin@# =
  log "Painting qualifiedAssociation skin...";

  fill objectBox(@#) withcolor @#info.background;

  draw objectBox(@#)
       withpen pencircle scaled @#info.borderScale withcolor @#info.borderColor;
enddef;

vardef QualifiedAssociation_draw@# =
  PictureStack_layout.@#elementsStack;

  QualifiedAssociation_paintSkin.@#;
  drawObject(@#elementsStack);
enddef;
