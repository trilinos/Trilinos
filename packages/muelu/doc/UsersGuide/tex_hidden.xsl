<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text"/>

  <xsl:template match="/masterlist">
    <xsl:for-each select="/masterlist/*/parameter">

      <!-- If default field is present, include it -->
      <xsl:choose>
        <xsl:when test="default">
\cbb{<xsl:value-of select="name"/>}{<xsl:value-of select="type"/>}{<xsl:value-of select="default"/>}{<xsl:value-of select="description"/>}
        </xsl:when><xsl:otherwise>
\cba{<xsl:value-of select="name"/>}{<xsl:value-of select="type"/>}{<xsl:value-of select="description"/>}
        </xsl:otherwise>
      </xsl:choose>

    </xsl:for-each>
  </xsl:template>

</xsl:stylesheet>
