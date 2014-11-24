<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text"/>

  <xsl:param name="section"/>

  <xsl:template match="/masterlist">
    <xsl:for-each select="/masterlist/*[local-name()=$section]/parameter">

      <!-- Skip parameter if visible field is not true. If it is absent, or is true, do the transformation -->
      <xsl:if test="compatibility-ML">

        <!-- If default field is present, include it -->
<xsl:choose>
        <xsl:when test="name-ML">
\mlcbb{<xsl:value-of select="name-ML"/>
        </xsl:when>
        <xsl:otherwise> 
\mlcbb{<xsl:value-of select="name"/>
      </xsl:otherwise>
</xsl:choose>}{<xsl:choose>
        <xsl:when test="type-ML">
  <xsl:value-of select="type-ML"/>
        </xsl:when>
        <xsl:otherwise> 
  <xsl:value-of select="type"/>
      </xsl:otherwise>
</xsl:choose>}{<xsl:choose>
      <xsl:when test="default-ML">
  <xsl:value-of select="default-ML"/>
      </xsl:when><xsl:otherwise>
  <xsl:value-of select="default"/>
      </xsl:otherwise>
</xsl:choose>}{<xsl:choose>
      <xsl:when test="compatibility-ML">
  <xsl:value-of select="compatibility-ML"/>
      </xsl:when>
</xsl:choose>}{<xsl:choose>
      <xsl:when test="description-ML">
  <xsl:value-of select="description-ML"/>
      </xsl:when><xsl:otherwise>    
  <xsl:value-of select="description"/>}
      </xsl:otherwise>
</xsl:choose>}   
      </xsl:if>

    </xsl:for-each>
  </xsl:template>

</xsl:stylesheet>
