<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" indent="no"/>

  <xsl:template match="/facade">

  <xsl:value-of select="/facade/interpreter-logic"/>

</xsl:template>

</xsl:stylesheet>