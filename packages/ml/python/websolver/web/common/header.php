<?php 

  $styleSheetString = "";
  $altStyleSheetString = "";
  $filename = "";
  $title = "";
  $styles = "";
  $scripts = "";
  $bodyAttributes = "";
  $dir = "";

  function setFilename ($in_filename) { global $filename; $filename = $in_filename; }   
  function setTitle ($in_title) { global $title; $title = $in_title; }    
  function setStyles ($in_styles) { global $styles; $styles = $in_styles; }    
  function setScripts ($in_scripts) { global $scripts; $scripts = $in_scripts; }  
  function setBodyAttributes ($in_bodyAttributes) { global $bodyAttributes; $bodyAttributes = $in_bodyAttributes; }  
  function setDir ($in_dir) { global $dir; $dir = $in_dir; } 
  
  function includeStyleSheet ($in_styleSheet) { 
    global $styleSheetString; 
    $styleSheetString .= "<link rel=\"stylesheet\" type=\"text/css\" href=\"$in_styleSheet\" \\>\n"; 
  }
  
  function includeAltStyleSheet ($in_altStyleSheet, $in_title) { 
    global $altStyleSheetString; 
    $altStyleSheetString .= "<link rel=\"alternate stylesheet\" type=\"text/css\" href=\"$in_altStyleSheet\" title=\"$in_title\" \\>\n"; 
  }
  
?>