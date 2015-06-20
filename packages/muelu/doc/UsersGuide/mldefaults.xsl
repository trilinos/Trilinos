<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" indent="no"/>

  <xsl:template match="/masterlist">"&lt;ParameterList name=\"MueLu\"&gt;"<xsl:for-each select="/masterlist/*/parameter">
  
  <xsl:choose>
      
      <!-- automatically process ML parameters only if there is no comment-ML parameter!
           Parameters with a comment-ML are supposed to need some hard-coded special handling
           in the source code! -->
      <xsl:when test="comment-ML">
      </xsl:when> 
      <xsl:otherwise>&lt;Parameter name=\"<xsl:choose>
        <xsl:when test="name-ML">
  <xsl:value-of select="name-ML"/>
        </xsl:when>
        <xsl:otherwise> 
  <xsl:value-of select="name"/>
      </xsl:otherwise>
      </xsl:choose>\" type=\"<xsl:choose>
      <xsl:when test="type-ML">
<xsl:value-of select="type-ML"/>
      </xsl:when>
      <xsl:otherwise>
<xsl:value-of select="type"/>
      </xsl:otherwise>
      </xsl:choose>\" value=\"<xsl:choose>
      <xsl:when test="default-ML">
<xsl:value-of select="default-ML"/>
      </xsl:when>
      <xsl:otherwise>
<xsl:value-of select="default"/>
      </xsl:otherwise>
      </xsl:choose>\"&gt;
      </xsl:otherwise>
    </xsl:choose>
    </xsl:for-each>
"&lt;/ParameterList&gt;"
</xsl:template>

</xsl:stylesheet>