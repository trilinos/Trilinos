<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" indent="no"/>

  <xsl:template match="/masterlist">

  <!-- Automatically generate xml generator C++ fragments -->
  <!-- The code maps the MueLu name to the corresponding XML string with a variable as value.
       The XML format corresponds to the simple XML input deck of MueLu -->
  <!-- The C++ code snippets are used for the ML 2 MueLu parameter translation class. Therefore we only process data
       which is not marked by an "comment-ML" entry. These parameters are supposed to need some hard-coded special handling
       for the parameter transition from ML to MueLu! -->
  <xsl:for-each select="/masterlist/*/parameter">
    <xsl:choose>
      <xsl:when test="comment-ML">
      </xsl:when>
      <xsl:otherwise>
    if (name == "<xsl:value-of select="name"/>") { ss &lt;&lt; "&lt;Parameter name=\"<xsl:value-of select="name"/>\" type=\"<xsl:value-of select="type"/>\" value=" &lt;&lt; value &lt;&lt; "/&gt;"; return ss.str(); }      </xsl:otherwise>
      </xsl:choose>
  </xsl:for-each>

</xsl:template>

</xsl:stylesheet>