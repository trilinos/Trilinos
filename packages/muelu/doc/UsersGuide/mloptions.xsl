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
          <xsl:when test="default-ML">
\mlcbb{<xsl:value-of select="name-ML"/>}{<xsl:value-of select="type-ML"/>}{<xsl:value-of select="default-ML"/>}{<xsl:value-of select="compatibility-ML"/>}{<xsl:value-of select="description-ML"/>}
          </xsl:when><xsl:otherwise>
\mlcba{<xsl:value-of select="name-ML"/>}{<xsl:value-of select="type-ML"/>}{<xsl:value-of select="compatibility-ML"/>}{<xsl:value-of select="description-ML"/>}
          </xsl:otherwise>
        </xsl:choose>

      </xsl:if>

    </xsl:for-each>
  </xsl:template>

</xsl:stylesheet>
