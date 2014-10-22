<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:param name="prob_type"/>
  <xsl:output method="text" indent="no"/>

  <xsl:template match="/masterlist">
    "&lt;ParameterList name=\"MueLu\"&gt;"
    <xsl:for-each select="/masterlist/*/parameter">
      <xsl:if test="*[local-name()=$prob_type]">
        <xsl:choose>
          <!-- handles "prob_type" parameter sublists -->
          <xsl:when test="type='\parameterlist'">
            <xsl:if test="./*[name() = $prob_type]">
              <xsl:call-template name="ParseProbSpecificSubList"/>
            </xsl:if>
          </xsl:when>
          <!-- handles "prob_type" params -->
          <xsl:otherwise>
            "&lt;Parameter name=\"<xsl:value-of select="name"/>\" type=\"<xsl:value-of select="type"/>\" value=\"<xsl:value-of select="*[local-name()=$prob_type]"/>\"/&gt;"
          </xsl:otherwise>
        </xsl:choose>
      </xsl:if>
    </xsl:for-each>
    "&lt;/ParameterList&gt;"
  </xsl:template>

  <xsl:template name="ParseProbSpecificSubList">
    "&lt;ParameterList name=\"<xsl:value-of select="name"/>\"&gt;"
    <xsl:for-each select="./*[name() = $prob_type]/parameter">
      <xsl:choose>
        <!-- the problem-specific sublist can contain a sublist -->
        <xsl:when test="type='\parameterlist'">
          <xsl:call-template name="ParseSubList"/>
        </xsl:when>
        <xsl:otherwise>
        "&lt;Parameter name=\"<xsl:value-of select="name"/>\" type=\"<xsl:value-of select="type"/>\" value=\"<xsl:value-of select="default"/>\"/&gt;"
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
    "&lt;/ParameterList&gt;"
  </xsl:template>

  <xsl:template name="ParseSubList">
    "&lt;ParameterList name=\"<xsl:value-of select="name"/>\"&gt;"
    <xsl:for-each select="./parameter">
      <xsl:choose>
        <!-- a sublist can contain a sublist -->
        <xsl:when test="type='\parameterlist'">
          <xsl:call-template name="ParseSubList"/>
        </xsl:when>
        <xsl:otherwise>
        "&lt;Parameter name=\"<xsl:value-of select="name"/>\" type=\"<xsl:value-of select="type"/>\" value=\"<xsl:value-of select="default"/>\"/&gt;"
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
    "&lt;/ParameterList&gt;"
  </xsl:template>



</xsl:stylesheet>
