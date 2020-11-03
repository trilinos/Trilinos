/*
  Copyright(C) 1999-2020 National Technology & Engineering Solutions
  of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
  NTESS, the U.S. Government retains certain rights in this software.

  See packages/seacas/LICENSE for details

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

      * Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.

      * Neither the name of NTESS nor the names of its
        contributors may be used to endorse or promote products derived
        from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
%{
#include <stdio.h>

#define YYPARSE_PARAM

#include "catalystParser.tab.hpp"

void track_location(YYLTYPE *loc, char *txt)
{
  loc->first_line = loc->last_line;
  loc->first_column = loc->last_column;
  int i;
  for(i = 0; txt[i] != '\0'; i++)
    {
    if(txt[i] == '\n')
      {
      loc->last_line++;
      loc->last_column = 0;
      }
    else
      {
      loc->last_column++;
      }
    }
}

#define YY_USER_ACTION track_location(yylloc, yytext);

%}

%option reentrant bison-bridge bison-locations
%option prefix="catalyst"

%%

 /* On Off Criteria Block */

(?i:onoffcriteria)                                                   { return(TOKONOFFCRITERIABLOCK); }
(?i:triggerthresholds)                                               { yylval->string=strdup(yytext); return(TOKTRIGGERTHRESHOLDS); }
(?i:experimental)                                                    { yylval->string=strdup(yytext); return(TOKEXPERIMENTAL); }
(?i:experimental\ +type)                                             { return(TOKEXPERIMENTALTYPE); }
(?i:threshold\ +variable\ +function)                                 { return(TOKTHRESHOLDVARIABLEFUNCTION); }
(?i:minimum)                                                         { yylval->string=strdup(yytext); return(TOKMINIMUM); }
(?i:maximum)                                                         { yylval->string=strdup(yytext); return(TOKMAXIMUM); }
(?i:mean)                                                            { yylval->string=strdup(yytext); return(TOKMEAN); }
(?i:sum)                                                             { yylval->string=strdup(yytext); return(TOKSUM); }
(?i:count)                                                           { yylval->string=strdup(yytext); return(TOKCOUNT); }
(?i:threshold\ +values)                                              { return(TOKTHRESHOLDVALUES); }
(?i:frames\ +after\ +trigger)                                        { return(TOKFRAMESAFTERTRIGGER); }
(?i:maximum\ +triggers\ +per\ +threshold)                           { return(TOKMAXIMUMTRIGGERSPERTHRESHOLD); }

 /* On Off Criteria Block */

 /* Experimental Block */

(?i:experimental\ +block)                                            { return(TOKEXPERIMENTALBLOCK); }

 /* Experimental Block */

 /* Visual Marker */

(?i:visualmarker)                                                    { return(TOKVISUALMARKER); }
(?i:shape)                                                           { return(TOKSHAPE); }
(?i:color)                                                           { return(TOKCOLOR); }
(?i:absolute\ +size)                                                 { return(TOKABSOLUTESIZE); }
(?i:relative\ +size)                                                 { return(TOKRELATIVESIZE); }
(?i:orientation)                                                     { return(TOKORIENTATION); }
(?i:resolution)                                                      { return(TOKRESOLUTION); }
(?i:marker\ +at\ +relative\ +point)                                  { return(TOKMARKERATRELATIVEPOINT); }
(?i:marker\ +at\ +absolute\ +point)                                  { return(TOKMARKERATABSOLUTEPOINT); }
(?i:marker\ +at\ +node)                                              { return(TOKMARKERATNODE); }
(?i:marker\ +at\ +element)                                           { return(TOKMARKERATELEMENT); }
(?i:marker\ +at\ +data\ +point)                                      { return(TOKMARKERATDATAPOINT); }
(?i:vector\ +parameter)                                              { return(TOKVECTORPARAMETER); }
(?i:scalar\ +parameter)                                              { return(TOKSCALARPARAMETER); }
(?i:integer\ +parameter)                                             { return(TOKINTEGERPARAMETER); }
(?i:string\ +parameter)                                              { return(TOKSTRINGPARAMETER); }
(?i:min)                                                             { yylval->string=strdup(yytext); return(TOKMIN); }
(?i:max)                                                             { yylval->string=strdup(yytext); return(TOKMAX); }
(?i:none)                                                            { yylval->string=strdup(yytext); return(TOKNONE); }

 /* Visual Marker */

 /* Experimental Operation */

(?i:experimental\ +operation)                                        { return(TOKEXPERIMENTALOPERATION); }
(?i:type)                                                            { return(TOKTYPE); }
(?i:vector\ +opparam)                                                { return(TOKVECTOROPPARAM); }
(?i:scalar\ +opparam)                                                { return(TOKSCALAROPPARAM); }
(?i:integer\ +opparam)                                               { return(TOKINTEGEROPPARAM); }
(?i:string\ +opparam)                                                { return(TOKSTRINGOPPARAM); }

 /* Experimental Operation */

 /* Calculator */

(?i:calculator)                                                      { return(TOKCALCULATOR); }
(?i:output\ +variable\ +name)                                        { return(TOKOUTPUTVARIABLENAME); }
(?i:function)                                                        { return(TOKFUNCTION); }
(?i:element\ +or\ +node\ +data)                                      { return(TOKELEMENTORNODEDATA); }
(?i:coordinate\ +result\ +flag)                                      { return(TOKCOORDINATERESULTFLAG); }
(?i:result\ +normals\ +flag)                                         { return(TOKRESULTNORMALSFLAG); }
(?i:python\ +calculator\ +flag)                                      { return(TOKPYTHONCALCULATORFLAG); }

 /* Calculator */

 /* Warp By Vector */

(?i:warpbyvector)                                                    { return(TOKWARPBYVECTOR); }
(?i:variable\ +vector)                                               { return(TOKVARIABLEVECTOR); }

 /* Warp By Vector */

 /* Transform */

(?i:transform)                                                       { return(TOKTRANSFORM); }
(?i:rotate)                                                          { return(TOKROTATE); }
(?i:scale)                                                           { return(TOKSCALE); }
(?i:translate)                                                       { return(TOKTRANSLATE); }

 /* Transform */

 /* Extract Block */

(?i:extractblock)                                                    { return(TOKEXTRACTBLOCK); }
(?i:include\ +blocks)                                                { return(TOKINCLUDEBLOCKS); }
(?i:exclude\ +blocks)                                                { return(TOKEXCLUDEBLOCKS); }

 /* Extract Block */

 /* Group */

(?i:operation\ +group\ +list)                                        { return(TOKOPERATIONGROUPLIST); }
(?i:group)                                                           { return(TOKGROUP); }

 /* Group */

 /* Plot Over Time */;

(?i:plot\ +id)                                                       { return(TOKPLOTID); }
(?i:plot\ +mean)                                                     { return(TOKPLOTMEAN); }
(?i:plot\ +maximum)                                                  { return(TOKPLOTMAXIMUM); }
(?i:plot\ +minimum)                                                  { return(TOKPLOTMINIMUM); }
(?i:plot\ +over\ +time)                                              { return(TOKPLOTOVERTIME); }

 /* Plot Over Time */;

 /* Scatter Plot */;

(?i:plot\ +basedirectory)                                            { return(TOKPLOTBASEDIRECTORY); }
(?i:plot\ +basename)                                                 { return(TOKPLOTBASENAME); }
(?i:y\ +axis\ +maximum\ +range)                                      { return(TOKYAXISMAXIMUMRANGE); }
(?i:x\ +axis\ +maximum\ +range)                                      { return(TOKXAXISMAXIMUMRANGE); }
(?i:y\ +axis\ +minimum\ +range)                                      { return(TOKYAXISMINIMUMRANGE); }
(?i:x\ +axis\ +minimum\ +range)                                      { return(TOKXAXISMINIMUMRANGE); }
(?i:y\ +axis\ +use\ +cumulative\ +data\ +range)                      { return(TOKYAXISUSECUMULATIVEDATARANGE); }
(?i:x\ +axis\ +use\ +cumulative\ +data\ +range)                      { return(TOKXAXISUSECUMULATIVEDATARANGE); }
(?i:y\ +axis\ +use\ +current\ +data\ +range)                         { return(TOKYAXISUSECURRENTDATARANGE); }
(?i:x\ +axis\ +use\ +current\ +data\ +range)                         { return(TOKXAXISUSECURRENTDATARANGE); }
(?i:y\ +axis\ +range)                                                { return(TOKYAXISRANGE); }
(?i:x\ +axis\ +range)                                                { return(TOKXAXISRANGE); }
(?i:y\ +axis\ +variable\ +tensor\ +component)                        { return(TOKYAXISVARIABLETENSORCOMPONENT); }
(?i:x\ +axis\ +variable\ +tensor\ +component)                        { return(TOKXAXISVARIABLETENSORCOMPONENT); }
(?i:y\ +axis\ +variable\ +vector\ +component)                        { return(TOKYAXISVARIABLEVECTORCOMPONENT); }
(?i:x\ +axis\ +variable\ +vector\ +component)                        { return(TOKXAXISVARIABLEVECTORCOMPONENT); }
(?i:y\ +axis\ +variable\ +vector\ +magnitude)                        { return(TOKYAXISVARIABLEVECTORMAGNITUDE); }
(?i:x\ +axis\ +variable\ +vector\ +magnitude)                        { return(TOKXAXISVARIABLEVECTORMAGNITUDE); }
(?i:y\ +axis\ +variable\ +scalar)                                    { return(TOKYAXISVARIABLESCALAR); }
(?i:x\ +axis\ +variable\ +scalar)                                    { return(TOKXAXISVARIABLESCALAR); }
(?i:all\ +variables)                                                 { return(TOKALLVARIABLES); }
(?i:scatter\ +plot)                                                  { return(TOKSCATTERPLOT); }

 /* Scatter Plot */;

 /* Contour */;

(?i:contour\ +value\ +sequence)                                      { return(TOKCONTOURVALUESEQUENCE); }
(?i:contour\ +value)                                                 { return(TOKCONTOURVALUE); }
(?i:end\ +contour)                                                   { return(TOKENDCONTOUR); }

 /* Contour */;

 /* Threshold */;

(?i:variable\ +tensor\ +component)                                   { return(TOKVARIABLETENSORCOMPONENT); }
(?i:variable\ +vector\ +component)                                   { return(TOKVARIABLEVECTORCOMPONENT); }
(?i:variable\ +vector\ +magnitude)                                   { return(TOKVARIABLEVECTORMAGNITUDE); }
(?i:variable\ +scalar)                                               { return(TOKVARIABLESCALAR); }
(?i:end\ +threshold)                                                 { return(TOKENDTHRESHOLD); }

 /* Threshold */;

 /* BloxClip */;

(?i:keep\ +inside\ +box)                                             { return(TOKKEEPINSIDEBOX); }
(?i:rotations)                                                       { return(TOKROTATIONS); }
(?i:relative\ +extents)                                              { return(TOKRELATIVEEXTENTS); }
(?i:absolute\ +extents)                                              { return(TOKABSOLUTEEXTENTS); }
(?i:center\ +at\ +element)                                           { return(TOKCENTERATELEMENT); }
(?i:center\ +at\ +node)                                              { return(TOKCENTERATNODE); }
(?i:center\ +at\ +relative\ +point)                                  { return(TOKCENTERATRELATIVEPOINT); }
(?i:center\ +at\ +absolute\ +point)                                  { return(TOKCENTERATABSOLUTEPOINT); }
(?i:end\ +boxclip)                                                   { return(TOKENDBOXCLIP); }
(?i:center\ +at\ data\ +point)                                       { return(TOKCENTERATDATAPOINT); }

 /* BoxClip */;

 /* Slice */;

(?i:end\ +slice)                                                     { return(TOKENDSLICE); }

 /* Slice */;

 /* Clip */;

(?i:negative)                                                        { yylval->string=strdup(yytext); return(TOKNEGATIVE); }
(?i:positive)                                                        { yylval->string=strdup(yytext); return(TOKPOSITIVE); }
(?i:smooth)                                                          { yylval->string=strdup(yytext); return(TOKSMOOTH); }
(?i:crinkle)                                                         { yylval->string=strdup(yytext); return(TOKCRINKLE); }
(?i:side\ +to\ +keep)                                                { return(TOKSIDETOKEEP); }
(?i:cut\ +type)                                                      { return(TOKCUTTYPE); }
(?i:plane\ +normal)                                                  { return(TOKPLANENORMAL); }
(?i:element\ +on\ +plane)                                            { return(TOKELEMENTONPLANE); }
(?i:node\ +on\ +plane)                                               { return(TOKNODEONPLANE); }
(?i:relative\ +point\ +on\ +plane)                                   { return(TOKRELATIVEPOINTONPLANE); }
(?i:absolute\ +point\ +on\ +plane)                                   { return(TOKABSOLUTEPOINTONPLANE); }
(?i:input)                                                           { return(TOKINPUT); }
(?i:end\ +clip)                                                      { return(TOKENDCLIP); }
(?i:data\ +point\ on\ +plane)                                        { return(TOKDATAPOINTONPLANE); }
(?i:plane\ +specification)                                           { return(TOKPLANESPECIFICATION); }
(?i:three\ +points)                                                  { return(TOKTHREEPOINTS); }
(?i:point\ +and\ +normal)                                            { return(TOKPOINTANDNORMAL); }
(?i:a)                                                               { yylval->string=strdup(yytext); return(TOKA); }
(?i:b)                                                               { yylval->string=strdup(yytext); return(TOKB); }
(?i:c)                                                               { yylval->string=strdup(yytext); return(TOKC); }

 /* Clip */;

 /* ImageSet */;

(?i:center)                                                          { yylval->string=strdup(yytext); return(TOKCENTER); }
(?i:point)                                                           { return(TOKPOINT); }
(?i:keep\ +outside)                                                  { yylval->string=strdup(yytext); return(TOKKEEPOUTSIDE); }
(?i:keep\ +inside)                                                   { yylval->string=strdup(yytext); return(TOKKEEPINSIDE); }
(?i:slice)                                                           { return(TOKSLICE); }
(?i:boxclip)                                                         { return(TOKBOXCLIP); }
(?i:extents)                                                         { return(TOKEXTENTS); }
(?i:normal)                                                          { return(TOKNORMAL); }
(?i:clip)                                                            { return(TOKCLIP); }
(?i:value\ +sequence)                                                { yylval->string=strdup(yytext); return(TOKVALUESEQUENCE); }
(?i:value\ +list)                                                    { yylval->string=strdup(yytext); return(TOKVALUELIST); }
(?i:contour)                                                         { return(TOKCONTOUR); }
(?i:keep\ +below)                                                    { yylval->string=strdup(yytext); return(TOKKEEPBELOW); }
(?i:keep\ +above)                                                    { yylval->string=strdup(yytext); return(TOKKEEPABOVE); }
(?i:keep\ +between)                                                  { yylval->string=strdup(yytext); return(TOKKEEPBETWEEN); }
(?i:tensor\ +component)                                              { yylval->string=strdup(yytext); return(TOKTENSORCOMPONENT); }
(?i:vector\ +component)                                              { yylval->string=strdup(yytext); return(TOKVECTORCOMPONENT); }
(?i:vector\ +magnitude)                                              { yylval->string=strdup(yytext); return(TOKVECTORMAGNITUDE); }
(?i:scalar)                                                          { yylval->string=strdup(yytext); return(TOKSCALAR); }
(?i:threshold)                                                       { return(TOKTHRESHOLD); }
(?i:axes\ +color)                                                    { return(TOKAXESCOLOR); }
(?i:text\ +color)                                                    { return(TOKTEXTCOLOR); }
(?i:edge\ +color)                                                    { return(TOKEDGECOLOR); }
(?i:surface\ +color)                                                 { return(TOKSURFACECOLOR); }
(?i:background\ +color)                                              { return(TOKBACKGROUNDCOLOR); }
(?i:image\ +size)                                                    { return(TOKIMAGESIZE); }
(?i:image\ +digit\ +count)                                           { return(TOKIMAGEDIGITCOUNT); }
(?i:png)                                                             { return(TOKPNG); }
(?i:jpg)                                                             { return(TOKJPG); }
(?i:image\ +format)                                                  { return(TOKIMAGEFORMAT); }
(?i:image\ +basedirectory)                                           { return(TOKIMAGEBASEDIRECTORY); }
(?i:image\ +basename)                                                { return(TOKIMAGEBASENAME); }
(?i:operation)                                                       { return(TOKOPERATION); }
(?i:operation2)                                                      { return(TOKOPERATION2); }
(?i:operation3)                                                      { return(TOKOPERATION3); }
(?i:imageset)                                                        { return(TOKIMAGESET); }
(?i:onoff\ +criteria)                                                { return(TOKONOFFCRITERIA); }
(?i:markers)                                                         { return(TOKMARKERS); }

 /* ImageSet */;

 /* Camera */;

(?i:camera\ +at\ +data\ +point\ +displaced)                          { return(TOKCAMERAATDATAPOINTDISPLACED); }
(?i:camera\ +at\ +data\ +point)                                      { return(TOKCAMERAATDATAPOINT); }
(?i:camera\ +at\ +node\ +displaced)                                  { return(TOKCAMERAATNODEDISPLACED); }
(?i:camera\ +at\ +element\ +displaced)                               { return(TOKCAMERAATELEMENTDISPLACED); }
(?i:camera\ +at\ +node)                                              { return(TOKCAMERAATNODE); }
(?i:camera\ +at\ +element)                                           { return(TOKCAMERAATELEMENT); }
(?i:camera\ +at\ +relative\ +point)                                  { return(TOKCAMERAATRELATIVEPOINT); }
(?i:camera\ +at\ +absolute\ +point)                                  { return(TOKCAMERAATABSOLUTEPOINT); }
(?i:look\ +direction)                                                { return(TOKLOOKDIRECTION); }
(?i:camera)                                                          { return(TOKCAMERA); }

 /* Camera */;

 /* MultiCamera8 */;

(?i:camera\ +fov)                                                    { return(TOKCAMERAFOV); }
(?i:up\ +vector)                                                     { return(TOKUPVECTOR); }
(?i:look\ +at\ +relative\ +distance)                                 { return(TOKLOOKATRELATIVEDISTANCE); }
(?i:look\ +at\ +absolute\ +distance)                                 { return(TOKLOOKATABSOLUTEDISTANCE); }
(?i:look\ +at\ +node)                                                { return(TOKLOOKATNODE); }
(?i:look\ +at\ +element)                                             { return(TOKLOOKATELEMENT); }
(?i:look\ +at\ +relative\ +point)                                    { return(TOKLOOKATRELATIVEPOINT); }
(?i:look\ +at\ +absolute\ +point)                                    { return(TOKLOOKATABSOLUTEPOINT); }
(?i:multicamera8)                                                    { return(TOKMULTICAMERA8); }

 /* MultiCamera8 */;

 /* Representation */;

(?i:highlight\ +subrange)                                            { return(TOKHIGHLIGHTSUBRANGE); }
(?i:color\ +legend\ +subrange)                                       { return(TOKCOLORLEGENDSUBRANGE); }
(?i:image\ +name\ +addon)                                            { return(TOKIMAGENAMEADDON); }
(?i:show\ +orientation\ +axes)                                       { return(TOKSHOWORIENTATIONAXES); }
(?i:show\ +z\ +axis\ +tic\ +mark\ +labels)                           { return(TOKSHOWZAXISTICMARKLABELS); }
(?i:show\ +y\ +axis\ +tic\ +mark\ +labels)                           { return(TOKSHOWYAXISTICMARKLABELS); }
(?i:show\ +x\ +axis\ +tic\ +mark\ +labels)                           { return(TOKSHOWXAXISTICMARKLABELS); }
(?i:show\ +z\ +axis\ +minor\ +tic\ +marks)                           { return(TOKSHOWZAXISMINORTICMARKS); }
(?i:show\ +y\ +axis\ +minor\ +tic\ +marks)                           { return(TOKSHOWYAXISMINORTICMARKS); }
(?i:show\ +x\ +axis\ +minor\ +tic\ +marks)                           { return(TOKSHOWXAXISMINORTICMARKS); }
(?i:show\ +z\ +axis\ +tic\ +marks)                                   { return(TOKSHOWZAXISTICMARKS); }
(?i:show\ +y\ +axis\ +tic\ +marks)                                   { return(TOKSHOWYAXISTICMARKS); }
(?i:show\ +x\ +axis\ +tic\ +marks)                                   { return(TOKSHOWXAXISTICMARKS); }
(?i:z\ +axis\ +label\ +name)                                         { return(TOKZAXISLABELNAME); }
(?i:y\ +axis\ +label\ +name)                                         { return(TOKYAXISLABELNAME); }
(?i:x\ +axis\ +label\ +name)                                         { return(TOKXAXISLABELNAME); }
(?i:show\ +z\ +axis\ +label)                                         { return(TOKSHOWZAXISLABEL); }
(?i:show\ +y\ +axis\ +label)                                         { return(TOKSHOWYAXISLABEL); }
(?i:show\ +x\ +axis\ +label)                                         { return(TOKSHOWXAXISLABEL); }
(?i:show\ +axes)                                                     { return(TOKSHOWAXES); }
(?i:time\ +annotation\ +position)                                    { return(TOKTIMEANNOTATIONPOSITION); }
(?i:show\ +time\ +annotation)                                        { return(TOKSHOWTIMEANNOTATION); }
(?i:invert\ +color\ +scale)                                          { return(TOKINVERTCOLORSCALE); }
(?i:preset\ +color\ +scale)                                          { return(TOKPRESETCOLORSCALE); }
(?i:color\ +legend\ +position)                                       { return(TOKCOLORLEGENDPOSITION); }
(?i:color\ +legend\ +maximum\ +range)                                { return(TOKCOLORLEGENDMAXIMUMRANGE); }
(?i:color\ +legend\ +minimum\ +range)                                { return(TOKCOLORLEGENDMINIMUMRANGE); }
(?i:color\ +legend\ +use\ +cumulative\ +data\ +range)                { return(TOKCOLORLEGENDUSECUMULATIVEDATARANGE); }
(?i:color\ +legend\ +use\ +current\ +data\ +range)                   { return(TOKCOLORLEGENDUSECURRENTDATARANGE); }
(?i:color\ +legend\ +range)                                          { return(TOKCOLORLEGENDRANGE); }
(?i:show\ +color\ +legend)                                           { return(TOKSHOWCOLORLEGEND); }
(?i:color\ +by\ +blockid)                                            { return(TOKCOLORBYBLOCKID); }
(?i:color\ +by\ +solid\ +color)                                      { return(TOKCOLORBYSOLIDCOLOR); }
(?i:color\ +by\ +tensor\ +component)                                 { return(TOKCOLORBYTENSORCOMPONENT); }
(?i:color\ +by\ +vector\ +component)                                 { return(TOKCOLORBYVECTORCOMPONENT); }
(?i:color\ +by\ +vector\ +magnitude)                                 { return(TOKCOLORBYVECTORMAGNITUDE); }
(?i:color\ +by\ +scalar)                                             { return(TOKCOLORBYSCALAR); }
(?i:show\ +bounding\ +box)                                           { return(TOKSHOWBOUNDINGBOX); }
(?i:show\ +edges)                                                    { return(TOKSHOWEDGES); }
(?i:show\ +surfaces)                                                 { return(TOKSHOWSURFACES); }
(?i:top)                                                             { yylval->string=strdup(yytext); return(TOKTOP); }
(?i:left)                                                            { yylval->string=strdup(yytext); return(TOKLEFT); }
(?i:bottom)                                                          { yylval->string=strdup(yytext); return(TOKBOTTOM); }
(?i:right)                                                           { yylval->string=strdup(yytext); return(TOKRIGHT); }
(?i:top\ +left)                                                      { yylval->string=strdup(yytext); return(TOKTOPLEFT); }
(?i:top\ +right)                                                     { yylval->string=strdup(yytext); return(TOKTOPRIGHT); }
(?i:bottom\ +left)                                                   { yylval->string=strdup(yytext); return(TOKBOTTOMLEFT); }
(?i:bottom\ +right)                                                  { yylval->string=strdup(yytext); return(TOKBOTTOMRIGHT); }
(?i:representation)                                                  { return(TOKREPRESENTATION); }
(?i:representation2)                                                 { return(TOKREPRESENTATION2); }
(?i:representation3)                                                 { return(TOKREPRESENTATION3); }

 /* Representation */;

 /* Catalyst */;

(?i:catalyst)                                                        { return(TOKCATALYST); }
(?i:begin)                                                           { return(TOKBEGIN); }
(?i:end)                                                             { return(TOKEND); }

 /* Catalyst */;

"="                                                                  { return(TOKEQUAL); }
(?i:node)                                                            { yylval->string=strdup(yytext); return(TOKNODE); }
(?i:element)                                                         { yylval->string=strdup(yytext); return(TOKELEMENT); }
(?i:false)                                                           { yylval->integer=0; return(TOKFALSE); }
(?i:true)                                                            { yylval->integer=1; return(TOKTRUE); }
[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?                               { yylval->floating=atof(yytext); return(TOKVALUE); }
[_\-a-zA-Z][a-zA-Z0-9_\.-\[\]\*\(\)\?!&\^%@~\+\\]*                   { yylval->string=strdup(yytext); return(TOKIDENTIFIER); }
\"(\\.|[^\\"\n])*\"                                                  { yylval->string=strdup(yytext); return(TOKSTRING); }
[#|$].*                                                              { }
[ \t\v\f\n]+                                                         { }

%%

