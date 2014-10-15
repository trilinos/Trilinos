<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:param name="prob_type"/>
  <xsl:output method="text" indent="no"/>
  <xsl:template match="/masterlist">"&lt;ParameterList name=\"MueLu\"&gt;"
    <xsl:for-each select="/masterlist/*/parameter">
        <xsl:if test="*[local-name()=$prob_type]">
          "&lt;Parameter name=\"<xsl:value-of select="name"/>\" type=\"<xsl:value-of select="type"/>\" value=\"<xsl:value-of select="*[local-name()=$prob_type]"/>\"/&gt;"
        </xsl:if>
    </xsl:for-each>
"&lt;/ParameterList&gt;"
</xsl:template>
</xsl:stylesheet>
