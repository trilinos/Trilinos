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

if known _metauml_class_mp:
  expandafter endinput
fi;
_metauml_class_mp:=1;

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

string accessPublic, accessProtected, accessPrivate;
accessPublic:="+";
accessProtected:="#";
accessPrivate:="-";

vardef ClassInfo@#=
  attributes(@#);
  var(numeric) featureAccessSpacing, accessWidth, accessHeight;
  var(color)   foreColor, borderColor;

  FontInfo.@#nameFont(metauml_defaultFont, defaultscale);

  FontInfo.@#attributeFont (metauml_defaultFont, defaultscale);
  FontInfo.@#methodFont    (metauml_defaultFont, defaultscale);
  FontInfo.@#stereotypeFont(metauml_defaultFont, .7);

  ShadeInfo.@#iShade;

  @#featureAccessSpacing := 2;
  @#accessWidth := 7;
  @#accessHeight := 6.75;

  @#foreColor := .9white;
  @#borderColor := black;

  PictureInfo.@#iName      (2, 2, 1, 3)(@#nameFont);

  PictureInfo.@#iStereotype(2, -5, 2, 2)(@#stereotypeFont);
  PictureInfo.@#iAttribute (2, 2, 1.25, 0)(@#attributeFont);
  PictureInfo.@#iMethod    (2, 2, 1.25, 0)(@#methodFont);

  PictureStackInfo.@#iStereotypeStack(2, 2, 1, 1)(5.5)(@#iStereotype);
  PictureStackInfo.@#iAttributeStack (2, 2, 2.5, 2)(10.5)(@#iAttribute); % 8.5
  PictureStackInfo.@#iMethodStack    (2, 2, 2.5, 2)(10.5)(@#iMethod);
  @#iAttributeStack.iPict.bottom := 2;
  %@#iAttributeStack.iPict.boxed := 1;
  @#iMethodStack.iPict.bottom := 2;
  %@#iMethodStack.iPict.boxed := 1;

  LockInfo.@#iLock(@#accessWidth, @#accessHeight, .6, .15, .55, .4white, .6white, .7white, .3white);

  @#iName.ignoreNegativeBase := 1;
  @#iAttributeStack.iPict.ignoreNegativeBase := 1;
  @#iMethodStack.iPict.ignoreNegativeBase := 1;
  @#iStereotypeStack.iPict.ignoreNegativeBase := 1;
enddef;

vardef ClassInfoCopy@#(text src)=
  log "ClassInfoCopy: Copying class";

  attributes(@#);
  var(numeric) featureAccessSpacing, accessWidth, accessHeight;
  var(color)   foreColor, borderColor;

  FontInfoCopy.@#nameFont(src.nameFont);

  FontInfoCopy.@#attributeFont(src.attributeFont);
  FontInfoCopy.@#methodFont(src.methodFont);
  FontInfoCopy.@#stereotypeFont(src.stereotypeFont);

  ShadeInfoCopy.@#iShade(src.iShade);

  @#featureAccessSpacing := src.featureAccessSpacing;
  @#accessWidth := src.accessWidth;
  @#accessHeight := src.accessHeight;

  @#foreColor := src.foreColor;
  @#borderColor := src.borderColor;

  PictureInfoCopy.@#iName      (src.iName);

  PictureInfoCopy.@#iStereotype(src.iStereotype);
  PictureInfoCopy.@#iAttribute (src.iAttribute);
  PictureInfoCopy.@#iMethod    (src.iMethod);

  PictureStackInfoCopy.@#iStereotypeStack(src.iStereotypeStack);
  PictureStackInfoCopy.@#iAttributeStack (src.iAttributeStack);
  PictureStackInfoCopy.@#iMethodStack    (src.iMethodStack);

  LockInfoCopy.@#iLock(src.iLock);
enddef;

ClassInfo.iClass;

ClassInfoCopy.iClassNameOnly(iClass);
iClassNameOnly.iName.top    := 10;
iClassNameOnly.iName.bottom := 10;
iClassNameOnly.iAttributeStack.top := 0;
iClassNameOnly.iAttributeStack.bottom := 0;
iClassNameOnly.iMethodStack.top := 0;
iClassNameOnly.iMethodStack.bottom := 0;

ClassInfoCopy.iInterface(iClass);
iInterface.iAttributeStack.top := 0;
iInterface.iAttributeStack.bottom := 0;
iInterface.iName.iFont.name := metauml_defaultFontOblique;

ClassInfoCopy.iAbstractClass(iInterface);

%
% CLASS
%
vardef defClass@#(expr pname) =
  ObjectEquations(@#);
  @#className := "Class";

  string @#name;
  @#name = pname;

  numeric @#nStereotypes;
  string @#stereotypes[];

  string @#attributes[];
  string @#attributesAccess[];

  string @#methods[];
  string @#methodsAccess[];

  numeric @#nAttrs;
  numeric @#nMethods;

  @#nStereotypes := 0;
  @#nAttrs   := 0;
  @#nMethods := 0;
enddef;

vardef addAttribute@#(expr pcontent) =
  string access;
  access := substring(0,1) of pcontent;
  if (not (access = accessPublic))  and
     (not (access = accessPrivate)) and
     (not (access = accessProtected)):

    @#.attributes[@#.nAttrs] := pcontent;
    @#.attributesAccess[@#.nAttrs] := accessProtected;
  else:
    save from;
    from := 1;
    if (substring(1,2) of pcontent) = " ": from := 2; fi;

    @#.attributes[@#.nAttrs] := substring(from, infinity) of pcontent;
    @#.attributesAccess[@#.nAttrs] := access;
  fi;

  @#.nAttrs := @#.nAttrs + 1;
enddef;

vardef addMethod@#(expr pcontent) =
  string access;
  access := substring(0,1) of pcontent;
  if (not (access = accessPublic))  and
     (not (access = accessPrivate)) and
     (not (access = accessProtected)):

    @#.methods[@#.nMethods] := pcontent;
    @#.methodsAccess[@#.nMethods] := accessProtected;
  else:
    @#.methods[@#.nMethods] := substring(1, 999) of pcontent;
    @#.methodsAccess[@#.nMethods] := access;
  fi;

  @#.nMethods := @#.nMethods + 1;
enddef;

vardef classStereotype@#(expr pcontent) =
  @#stereotypes[@#nStereotypes] := pcontent;
  @#nStereotypes := @#nStereotypes + 1;
enddef;

vardef classStereotypes@#(text stereotypes)=
  for stereotype = stereotypes:
    classStereotype@#(stereotype);
  endfor;
enddef;

vardef EClass@#(text _info)(expr name)(text attributes)(text methods)=
  log "EClass begin: " & str @#;
  defClass@#(name);

  log "Copying class info";
  ClassInfoCopy.@#info(_info);

  for a=attributes:
    log "Adding attribute ", a;
    addAttribute@#(a);
  endfor;
  for m=methods:
    log "Adding method ", m;
    addMethod@#(m);
  endfor;
enddef;

vardef Class@#(expr name)(text attributes)(text methods)=
  EClass@#(iClass)(name)(attributes)(methods);
enddef;

vardef Interface@#(expr name)(text methods)=
  EClass@#(iInterface)(name)()(methods);
enddef;

vardef EInterface@#(text _info)(expr name)(text methods)=
  EClass@#(_info)(name)()(methods);
enddef;

vardef ClassName@#(expr name)=
  EClass@#(iClassNameOnly)(name)()();
enddef;

vardef EClassName@#(text _info)(expr name)=
  EClass@#(_info)(name)()();
enddef;

vardef AbstractClass@#(expr name)(text attributes)(text methods)=
  EClass@#(iAbstractClass)(name)(attributes)(methods);
enddef;

vardef EAbstractClass@#(text _info)(expr name)(text methods)=
  EClass@#(_info)(name)()(methods);
enddef;

vardef Class_border@#=
  objectBox(@#)
enddef;

vardef Class_layout@# =
  if @#layedout = 1:
    log "Class " & (str @#) & " has already been layed out";
  else:
    @#layedout := 1;
    log "Class layout: " & (str @#);

    @#maxFeatureWidth := 0;

    EPictureStack.@#stereotypeStack(@#info.iStereotypeStack)
       (scantokens listArray(@#stereotypes)(@#nStereotypes))("vcenterbase");

    EPicture.@#namePict(@#info.iName)(@#name);
    EPictureStack.@#attributeStack(@#info.iAttributeStack)
       (scantokens listArray(@#attributes)(@#nAttrs))("vleftbase");
    EPictureStack.@#methodStack(@#info.iMethodStack)
       (scantokens listArray(@#methods)(@#nMethods))("vleftbase");

    layoutObjects(@#stereotypeStack, @#namePict, @#attributeStack, @#methodStack);

    log "Computing maxFeatureWidth";
    log @#stereotypeStack.width;
    log @#namePict.width;
    log @#attributeStack.width;
    log @#methodStack.width;
    log "...";
    @#maxFeatureWidth := lmax(@#stereotypeStack.width, @#namePict.width,
                              @#attributeStack.width, @#methodStack.width);
    log "Done computing maxFeatureWidth";

    log "max feature width: " & decimal @#maxFeatureWidth;

    @#namePict.midx = @#midx;
    @#stereotypeStack.midx = @#namePict.midx;

    @#stereotypeStack.top = @#top;
    @#namePict.top = @#stereotypeStack.bottom;
    @#attributeStack.top = @#namePict.bottom;
    @#methodStack.top = @#attributeStack.bottom;

    @#attributeStack.left = @#methodStack.left =
        @#left + @#info.accessWidth + @#info.featureAccessSpacing;

    @#width = @#maxFeatureWidth + @#info.accessWidth + @#info.featureAccessSpacing;
    @#bottom = @#methodStack.bottom;

    log "Class layout done...";
  fi;
enddef;

vardef Class_paintAccess@#(expr _access, _se) =
  save lock;
  ELock.lock(@#info.iLock)(_access);
  lock.se = _se;
  drawObject(lock);
enddef;

vardef Class_drawFeatures@#=
  log "Drawing stereotypes";
  drawObject(@#stereotypeStack);

  log "Drawing name";
  drawObject(@#namePict);
  log "Drawing attribute stack";
  drawObject(@#attributeStack);

  log "Drawing method stack";
  drawObject(@#methodStack);

  log "********************************** Drawing access locks";
  for i = 0 upto @#nAttrs-1:
    log @#attributesAccess[i];
    Class_paintAccess@#(@#attributesAccess[i],
      @#attributeStack.pict[i].sw + (-@#info.featureAccessSpacing, @#info.iAttributeStack.iPict.bottom));
  endfor;

  for i = 0 upto @#nMethods-1:
    log @#methodsAccess[i];
    Class_paintAccess@#(@#methodsAccess[i],
      @#methodStack.pict[i].sw + (-@#info.featureAccessSpacing, @#info.iMethodStack.iPict.bottom));
  endfor;
enddef;

vardef Class_paintSkin@# =
  log "Painting class skin...";
  nameAttributeY := @#attributeStack.top;
  attributeMethodY := @#methodStack.top;

 %fill Class_border.@# shifted (@#info.shade,-@#info.shade)
 %     withcolor .7white withpen currentpen scaled 1;
 drawObjectShade(@#);

 fill Class_border.@# withcolor @#info.foreColor;
 draw Class_border.@# withcolor @#info.borderColor;

 draw (xpart @#nw, nameAttributeY)--(xpart @#se, nameAttributeY) withcolor @#info.borderColor;
 draw (xpart @#nw, attributeMethodY)--(xpart @#se, attributeMethodY) withcolor @#info.borderColor;
enddef;

vardef Class_draw@#=
  log "draw class begin...";
  Class_layout.@#;
  objectEnsurePositioning.@#;
  Class_paintSkin.@#;
  Class_drawFeatures.@#;
enddef;

vardef Class_setDebugMode@#=
  @#.info.iName.boxed := 1;
  @#.info.iStereotypeStack.boxed := 1;
  @#.info.iStereotypeStack.iPict.boxed := 1;
  @#.info.iAttributeStack.boxed := 1;
  @#.info.iAttributeStack.iPict.boxed := 1;
  @#.info.iMethodStack.boxed := 1;
  @#.info.iMethodStack.iPict.boxed := 1;
enddef;
