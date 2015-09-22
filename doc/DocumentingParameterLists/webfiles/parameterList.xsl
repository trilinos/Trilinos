<xsl:stylesheet version="2.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="/">
    <html>
      <head>
        <link rel="stylesheet" type="text/css"
 	  href="common/parameterList/bootstrap/css/bootstrap.css" />
        <link rel="stylesheet" type="text/css"
   	  href="common/parameterList/trilinos-bootstrap.css" />

        <script type="text/javascript"
         src="common/parameterList/js/jquery.js"></script>
        <script type="text/javascript"
         src="common/parameterList/js/iframeResizer.contentWindow.min.js"></script>
        <script type="text/javascript"
         src="common/parameterList/bootstrap/js/bootstrap.js"></script>
        <script type="text/javascript">
   $('.accordion-toggle').click(function(){
           $("i", this).toggleClass("t-icon-arrow-right t-icon-arrow-down");
         });
        </script>

      </head>
      <body>
        <div class="collapse-style">
          <xsl:apply-templates select="ParameterList"/>
        </div>
      </body>
    </html>
  </xsl:template>


  <xsl:template match="ParameterList">

    <xsl:variable name="paramlistid" select="@id" />
    <ul>
    <li>
      <div class="accordion" id="accordion{$paramlistid}">
        <div class="accordion-group">
          <div class="accordion-heading">
            <a class="accordion-toggle" data-toggle="collapse"
              data-parent="#accordion{$paramlistid}"
              href="#collapse{$paramlistid}"><i class="t-icon-arrow-right"></i>
                <xsl:choose>
                  <xsl:when test="@name='ANONYMOUS'">Valid Parameters</xsl:when>
                  <xsl:otherwise>
                    <xsl:call-template name="quoteStr">
                      <xsl:with-param name="inputStr" select="@name"/>
                    </xsl:call-template>
                  </xsl:otherwise>
                </xsl:choose>
              (ParameterList)
            </a>
          </div><!--accordion-heading-->
          <div id="collapse{$paramlistid}" class="accordion-body collapse">
            <div class="accordion-inner"><xsl:apply-templates/></div>
          </div><!--id-->
        </div><!--accordion-group-->
      </div><!--accordion-->
    </li>
    </ul>
  </xsl:template>

  <xsl:template match="Parameter">

    <xsl:variable name="paramvalidatorid" select="@validatorId" />
    <xsl:variable name="default" select="@isDefault" />
    <xsl:variable name="paramid" select="concat(@validatorId,@id)" />
    <ul>
      <li>
      <div class="accordion" id="accordion{$paramid}">
        <div class="accordion-group">
          <div class="accordion-heading">
            <a class="accordion-toggle" data-toggle="collapse"
              data-parent="#accordion{$paramid}"
              href="#collapse{$paramid}"><i class="t-icon-arrow-right"></i>
                <xsl:call-template name="quoteStr">
                  <xsl:with-param name="inputStr" select="@name"/>
                </xsl:call-template>
              ( <xsl:value-of select="@type"/>
              <xsl:text>, default = </xsl:text>
              <xsl:choose>
                <xsl:when test="@type='string'">
                  <xsl:call-template name="quoteStr">
                    <xsl:with-param name="inputStr" select="@value"/>
                  </xsl:call-template>
                </xsl:when>

                <xsl:when test="@type='double'">
                  <xsl:call-template name="printFloat">
                    <xsl:with-param name="x" select="@value" />
                  </xsl:call-template>
                </xsl:when>

                 <xsl:otherwise><xsl:value-of select="@value"/></xsl:otherwise>
              </xsl:choose>
              )
            </a>
             </div><!--accordion-heading-->
          <div id="collapse{$paramid}" class="accordion-body collapse">
            <div class="accordion-inner">
            <div class="emph-italic"><xsl:value-of select="@docString"/></div>

            <xsl:choose>

              <xsl:when test="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@type='EnhancedNumberValidator(int)'">
                <xsl:variable name="min" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@min"
                  as="xs:integer"/>
                <xsl:variable name="max" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@max"
                  as="xs:integer"/>
                <xsl:variable name="step" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@step"
                  as="xs:integer"/>
                <p><span class="emph-bold">Valid range:</span> [<xsl:value-of select="$min"/> : <xsl:value-of select="$max"/>]</p>
                <xsl:if test="$step!='1'">
                  <p>(step = <xsl:value-of select="$step"/>)</p>
                </xsl:if>
              </xsl:when>

              <xsl:when test="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@type='ArrayValidator(EnhancedNumberValidator(int), int)'">
                <xsl:variable name="min" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/Validator/@min"
                  as="xs:integer"/>
                <xsl:variable name="max" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/Validator/@max"
                  as="xs:integer"/>
                <xsl:variable name="step" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/Validator/@step"
                  as="xs:integer"/>
                <p><span class="emph-bold">Valid range:</span> [<xsl:value-of select="$min"/> : <xsl:value-of select="$max"/>]</p>
                <xsl:if test="$step!='1'">
                  <p>(step = <xsl:value-of select="$step"/>)</p>
                </xsl:if>
	      </xsl:when>

              <xsl:when test="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@type='EnhancedNumberValidator(double)'">
                <xsl:variable name="min" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@min"
                  as="xs:string"/>
                <xsl:variable name="max" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@max"
                  as="xs:string"/>
                <xsl:variable name="mymin">
                  <xsl:call-template name="printFloat">
                    <xsl:with-param name="x" select="$min"/>
                  </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="mymax">
                  <xsl:call-template name="printFloat">
                    <xsl:with-param name="x" select="$max"/>
                  </xsl:call-template>
                </xsl:variable>
                <p>[ <xsl:value-of select="$mymin"/> : <xsl:value-of select="$mymax"/> ]</p>
              </xsl:when>

              <xsl:when test="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@type='ArrayValidator(EnhancedNumberValidator(double), double)'">
                <xsl:variable name="min" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/Validator/@min"
                  as="xs:string"/>
                <xsl:variable name="max" select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/Validator/@max"
                  as="xs:string"/>
                <xsl:variable name="mymin">
                  <xsl:call-template name="printFloat">
                    <xsl:with-param name="x" select="$min"/>
                  </xsl:call-template>
                </xsl:variable>
                <xsl:variable name="mymax">
                  <xsl:call-template name="printFloat">
                    <xsl:with-param name="x" select="$max"/>
                  </xsl:call-template>
                </xsl:variable>
                <p>[ <xsl:value-of select="$mymin"/> : <xsl:value-of select="$mymax"/> ]</p>
              </xsl:when>

              <xsl:when test="contains(/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/@type,'StringIntegralValidator')">
                <table>
                  <tr>
                    <th>Valid strings</th>
                    <th>Description</th>
                  </tr>
                  <xsl:for-each select="/ParameterList/Validators/Validator[@validatorId=current()/@validatorId]/String">
                  <tr>
                    <td>
                      <xsl:call-template name="quoteStr">
                        <xsl:with-param name="inputStr" select="@stringValue"/>
                      </xsl:call-template>
                    </td>
                    <td>
                      <xsl:choose>
                        <xsl:when test="@stringDoc">
                          <xsl:value-of select="@stringDoc"/>
                        </xsl:when>
                        <xsl:otherwise>
                          <xsl:value-of select="@stringValue"/>
                        </xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </tr>
                </xsl:for-each>
                </table>
               </xsl:when>

              <xsl:otherwise>
              </xsl:otherwise>

            </xsl:choose>
            </div> <!-- accordion-inner -->
          </div> <!--id-->
        </div><!--accordion-group-->
      </div><!--accordion-->
    </li>
    </ul>

    <xsl:if test="ParameterList">
      <ul>
        <li><xsl:value-of select="@name"/></li>
        <xsl:apply-templates select="ParameterList" />
      </ul>
    </xsl:if>
  </xsl:template>

<!-- FUNCTION quoteStr -->
  <xsl:template name="quoteStr">
    <xsl:param name="inputStr" />
    <xsl:variable name="strInput" select="string($inputStr)" />
    <xsl:value-of select="concat('&#34;', $strInput,'&#34;')" />
  </xsl:template>
<!-- FUNCTION quoteStr -->

<!-- FUNCTION abs -->
  <xsl:template name="abs">
    <xsl:param name="x" />

    <xsl:choose>
      <xsl:when test="$x &lt; 0">
        <xsl:value-of select="0 - $x"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$x"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>
<!-- FUNCTION abs -->

<!-- FUNCTION power -->
  <xsl:template name="power">
    <xsl:param name="base" />
    <xsl:param name="pow" />

    <xsl:variable name="absPow">
      <xsl:call-template name="abs">
        <xsl:with-param name="x" select="$pow" />
      </xsl:call-template>
    </xsl:variable>

    <xsl:choose>

      <xsl:when test="$absPow = 0">
        <xsl:value-of select="1" />
      </xsl:when>

      <xsl:otherwise>

        <xsl:variable name="temp">
          <xsl:call-template name="power">
            <xsl:with-param name="base" select="$base" />
            <xsl:with-param name="pow" select="$absPow - 1" />
          </xsl:call-template>
        </xsl:variable>
        <xsl:choose>
          <xsl:when test="$pow = $absPow">
            <xsl:value-of select="$base * $temp" />
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="1 div ($base * $temp)" />
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>

    </xsl:choose>

  </xsl:template>
<!-- FUNCTION power -->

<!-- FUNCTION printFloat -->
 <xsl:template name="printFloat">
    <xsl:param name="x" />

    <xsl:variable name="mantissa" select="substring-before($x, 'e')" />
    <xsl:variable name="exp" select="substring-after($x, 'e')"  as="number" />
    <xsl:variable name="exponent">
      <xsl:choose>
        <xsl:when test="starts-with($exp, '+')">
          <xsl:value-of select="substring-after($exp, '+')" as="number" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="$exp" as="number" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="factor" select="1000" />
    <xsl:variable name="y" select="round($mantissa * $factor)" />
    <xsl:variable name="newMantissa" select="$y div $factor" />

<!--
    <p>mantissa: <xsl:value-of select="$mantissa"/></p>
    <p>exp: <xsl:value-of select="$exp"/></p>
    <p>exponent: <xsl:value-of select="$exponent"/></p>
    <p>factor: <xsl:value-of select="$factor"/></p>
    <p>y: <xsl:value-of select="$y"/></p>
    <p>newMantissa: <xsl:value-of select="$newMantissa"/></p>
-->
    <xsl:choose>

      <xsl:when test="$exponent &lt; -2 or $exponent &gt; 2" >
        <xsl:value-of select="$newMantissa" />e<xsl:value-of select="$exponent" />
      </xsl:when>

      <xsl:otherwise>
        <xsl:variable name="temp">
          <xsl:call-template name="power">
            <xsl:with-param name="base" select="10"/>
            <xsl:with-param name="pow" select="$exponent"/>
          </xsl:call-template>
        </xsl:variable>
       <xsl:value-of select="$newMantissa * $temp"/>
      </xsl:otherwise>

    </xsl:choose>
  </xsl:template>
<!-- FUNCTION printFloat -->

</xsl:stylesheet>
