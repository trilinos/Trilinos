// This JavaScript program is copyright by Jasper Bedaux (C) 2003.
// This script can be freely used and distributed as long as this
// copyright notice is preserved. If you use and or modify the script
// you have to mention the author in a reasonable visible way.
// THIS SCRIPT COMES WITHOUT ANY WARRANTY OF ANY KIND.
// SELLING OF THIS SCRIPT IS NOT ALLOWED.
// 
// A link to http://www.bedaux.net/cpp2html/ will be appreciated.

  keys = new Array() // containing C++ language elements
  keys.push({style:comment, start:/\s*\/\*[\s\S]*?\*\//mg})
  keys.push({style:comment, start:/\s*\/\//mg, end:/\n/mg, neglect:/\\|\?\?\//mg})
  keys.push({style:precompiler, start:/\s*?^\s*(?:#|\?\?=|%:)/mg, end:/\n/m, neglect:/\\[\s\S]|\?\?\/[\s\S]/m})
  keys.push({style:stringLiteral, start:/\s*(?:\bL)?"/mg, end:/"/m, neglect:/\\[\s\S]|\?\?\/[\s\S]/m})
  keys.push({style:charLiteral, start:/\s*(?:\bL)?'/mg, end:/'/m, neglect:/\\[\s\S]|\?\?\/[\s\S]/m})
  keys.push({style:floatLiteral, start:/\s*(?:(?:\b\d+\.\d*|\.\d+)(?:E[\+\-]?\d+)?|\b\d+E[\+\-]?\d+)[FL]?\b|\s*\b\d+\./mgi})
  keys.push({style:intLiteral, start:/\s*\b(?:0[0-7]*|[1-9]\d*|0x[\dA-F]+)(?:UL?|LU?)?\b/mgi})
  keys.push({style:boolLiteral, start:/\s*\b(?:true|false)\b/mg})
  keys.push({style:types, start:/\s*\b(?:bool|char|double|float|int|long|short|signed|unsigned|void|wchar_t)\b/mg})
  keys.push({style:flowControl, start:/\s*\b(?:break|case|catch|continue|default|do|else|for|goto|if|return|switch|throw|try|while)\b/mg})
  keys.push({style:keyword, start:/\s*\b(?:asm|auto|class|const_cast|const|delete|dynamic_cast|enum|explicit|export|extern|friend|inline|main|mutable|namespace|new|operator|private|protected|public|register|reinterpret_cast|sizeof|static_cast|static|struct|template|this|typedef|typeid|typename|union|using|virtual|volatile|and_eq|and|bitand|bitor|compl|not_eq|not|or_eq|or|xor_eq|xor)\b/mg})
  keys.push({style:operator, start:/\s*[\{\}\[\]\(\)<>%:;\.\?\*\+\-\^&\|~!=,\\]+|\s*\//mg})

  function setOptions(useXhtml) {
    windowSettings = "menubar,scrollbars,status,resizable,height=600,width=600"
    if (!useXhtml) {
      outputStart = "<HTML>\n<HEAD>\n<!" +
        "-- GENERATED WITH OPTION \"USE_OLD_HTML_INSTEAD_OF_VALID_XHTML_1.1\" --" +
        ">\n<TITLE>C++ code colored by C++2HTML</TITLE>\n" +
        "<META NAME=\"generator\" CONTENT=\"C++2HTML by Jasper Bedaux\">\n<" +
        "!-- To generate your own colored code visit http://www.bedaux.net/cpp2html/ --" +
        ">\n</HEAD>\n<BODY BGCOLOR=\"#FFFFFF\">\n<PRE>"
      outputEnd = "</PRE>\n</BODY>\n</HTML>\n"
      for (var i = 0; i != keys.length; ++i) { // set xhtml tags
        keys[i].before = ""
        if (keys[i].style.bold) keys[i].before += "<B>"
        if (keys[i].style.italic) keys[i].before += "<I>"
        keys[i].before += "<FONT COLOR=\"" + keys[i].style.color + "\">"
        keys[i].after = "</FONT>"
        if (keys[i].style.italic) keys[i].after += "</I>"
        if (keys[i].style.bold) keys[i].after += "</B>"
      }
    }
    else {
      outputStart = "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" \"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">\n" +
        "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\">\n" +
        "<head>\n<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />\n" +
        "<title>C++ code colored by C++2HTML</title>\n" +
        "<meta name=\"generator\" content=\"C++2HTML by Jasper Bedaux\" />\n<" +
        "!-- To generate your own colored code visit http://www.bedaux.net/cpp2html/ --" +
        ">\n<style type=\"text/css\">\n"
      for (var i = 1; i != keys.length; ++i) {
        outputStart += "." + keys[i].style.name + " { color: " + keys[i].style.color + "; "
        if (keys[i].style.bold) outputStart += "font-weight: bold; "
        if (keys[i].style.italic) outputStart += "font-style: italic; "
        outputStart += "}\n"
      }
      outputStart += "</style>\n</head>\n<body>\n<pre>"
      outputEnd = "</pre>\n</body>\n</html>\n"
      for (var i = 0; i != keys.length; ++i) { // set xhtml tags
        keys[i].before = "<span class=\"" + keys[i].style.name + "\">"
        keys[i].after = "</span>"
      }
    }
  }

  function toHTML(s) { // convert special chars
    s = s.split("&").join("&amp;");
    s = s.split("<").join("&lt;");
    return s.split(">").join("&gt;");
  }

function colorCode(f) {
  if (versionOK) {
    setOptions(true);
    var s = f.script.value // get code
    s = s.replace(/\r\n?/gm, "\n") // convert DOS and/or MAC to UNIX
    var out = window.open("", "", windowSettings)
    var keyString = ""
    var match = 0
    out.document.open()
    out.document.write(outputStart)

    var previousMatch = -1
    for (var i = 0; i != keys.length; ++i) keys[i].startPos = -1
    for (var position = 0; position != s.length; position = keys[match].endPos) {
      for (var i = 0; i != keys.length; ++i) {
        if (keys[i].startPos < position) { // update needed
          keys[i].start.lastIndex = position
          var result = keys[i].start.exec(s)
          if (result != null) {
            keys[i].startPos = result.index
            keys[i].endPos = keys[i].start.lastIndex
          }
          else keys[i].startPos = keys[i].endPos = s.length
        }
      }
      match = 0
      for (var i = 1; i < keys.length; ++i) // find first matching key
        if (keys[i].startPos < keys[match].startPos) match = i
      if (keys[match].end != undefined) { // end must be found
        var end = new RegExp(keys[match].end.source + "|" + keys[match].neglect.source, "mg")
        end.lastIndex = keys[match].endPos
        while (keys[match].endPos != s.length) {
          result = end.exec(s)
          if (result != null) {
            if (result[0].search(keys[match].end) == 0) {
              keys[match].endPos = end.lastIndex
              break
            }
          }
          else keys[match].endPos = s.length
        }
      }
      var before = s.substring(position, keys[match].startPos)
      keyString = s.substring(keys[match].startPos, keys[match].endPos)
      var output = ""
      if ((before == "") && (match == previousMatch))
        output += toHTML(keyString)
      else {
        if (previousMatch != -1) output += keys[previousMatch].after 
        output += toHTML(before)
        if (keyString != "") output += keys[match].before + toHTML(keyString)
      }
      previousMatch = match
      out.document.write(output)
    }
    if (keyString != "") out.document.write(keys[match].after)
    out.document.write(outputEnd)
    out.document.close()
  }
  else alert("Sorry, your browser is too old. Minimum required: Microsoft Internet " +
    "Explorer 5.5, Mozilla 1.0, Netscape 6 or other Mozilla-based browser.\n")
}
