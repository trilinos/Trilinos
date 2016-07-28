<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" indent="no"/>

  <xsl:template match="/facade">

  <!-- Automatically generate xml generator C++ fragments -->
  <xsl:for-each select="/facade/input-defaults/parameter">"&lt;Parameter name=\"<xsl:value-of select="name"/>\" type=\"<xsl:value-of select="type"/>\" value=\"<xsl:value-of select="default"/>\"/&gt;"
  </xsl:for-each>

</xsl:template>

</xsl:stylesheet>