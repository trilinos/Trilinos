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

%% Usage of skins
%%
%% Skins are configuration files that, in geneal, customize the settings of 
%% MetaUML *Info objects. You most likely want to include the skin file after 
%% all other MetaUML files have been included:
%%
%% input metauml
%% input metauml_skin_simple
%% 
%% % normal usage of MetaUML follows here
%%
%% More advanced skin files might customize fonts for specific elements 
%% (eg for Class and Note), and also colors (eg foreground, shade, line etc.)

if known _metauml_skin_simple_mp: 
  expandafter endinput 
fi;
_metauml_skin_simple_mp:=1;

iClass.foreColor := blue;
iClass.iName.iFont.name := "cmr12";
iClass.iAttributeStack.iPict.iFont.name := "cmr10";
iClass.iMethodStack.iPict.iFont.name := "cmr10";


