#!/Net/local/gnu/bin/perl
#
# Script to email web-form data for inclusion in Sandia LGPL database.

##############################################################################
# Error-handling routine; print web page with error message.
sub dienow {
  my($msg) = @_;
  print "Content-type:text/html\n\n";
  print "<html><head><title>Download Error</title></head>\n";
  print "<body>\n";
  print "<p><p><h2>\n";
  print "Input Error:  \n";
  print "$msg<p>\n";
  print "</h2>\n";
  print "</body></html>\n";
  exit;
}

##############################################################################
# Main script

# Read CGI input.
use CGI qw(:all);
require "cgi-lib.pl";
&ReadParse;

# FOR TESTING ONLY
# Hash table for testing the script without CGI input.
#%in=("name","Joe Schmoe","company","SNL","addr1","MS111","addr2","PO Box 5800","city","ABQ","state","NM","zip","87185","email","jschmoe\@cs.sandia.gov","user_list","Y", "country","US");
# END FOR TESTING ONLY

############################################################################
# Instructions for using this script for a new software package.
#
# For each package...
#  + Create a unique download page that is just a UNIX link 
#    to the master download page; this link will be checked as
#    the referring page to determine which package was requested.
#    (This approach lets us know which package is requested
#    without having to set cookies.)
#  + Add an "if" test for the referring page to identify the
#    appropriate software package.  Test $prevpage against the
#    UNIX link file created above.
#  + Within the "if" test, set appropriate variables.  See the
#    $pkg=Zoltan case for an example.
#
# Variables that must be defined:
#  $pkg:          The name of the package.
#  $pkghome:      The home page for your software package.
#  $downloadfile  The URL for the file to be downloaded.
#  $recipient     Email addresses for people who should receive notices
#                 when the software is downloaded.  

$prevpage = $ENV{'HTTP_REFERER'};
if ($prevpage eq "http://www.cs.sandia.gov/Zoltan/Zoltan_download.html") {
  $pkg = "Zoltan";
  $pkghome = "http://www.cs.sandia.gov/Zoltan/";
  $downloadfile = "ftp://ftp.cs.sandia.gov/pub/papers/kddevin/zoltan_distrib.tar.gz";
  $recipient = "k_dragon\@attglobal.net";
}
#### Add more software packages here. ####
else {
  dienow "The software package name could not be determined."
}
############################################################################

# Error checking for required fields.
if ($in{"name"} eq "") {
  dienow "Please enter your name."
}
if ($in{"city"} eq "") {
  dienow "Please enter your city."
}
if ($in{"state"} eq "") {
  dienow "Please enter your state."
}
if ($in{"country"} eq "") {
  dienow "Please enter your country."
}

# Set flag for inclusion on user's list.
if ($in{"user_list"} eq "Y") {
  $ul = $in{"user_list"};
  if ($in{"email"} eq "") {
    dienow "You must enter your email address to be added to the User's mailing list.";
  }
# May later add code that automatically adds the user to the appropriate list.
}
else {
  $ul = "N";
}


# Send email to $database_mgr and $recipient when the software is downloaded.
# This email will be used to create the database for SNL Tech. Transfer.

$database_mgr = "kddevin\@cs.sandia.gov";
open(MAIL, "|/usr/lib/sendmail -i $database_mgr $recipient") || dienow "Unable to send mail.";
print MAIL "$pkg|$ul|$in{email}|$in{name}|$in{company}|$in{addr1}|$in{addr2}|$in{city}|$in{state}|$in{zip}|$in{country}|$ENV{'REMOTE_ADDR'}\n\n";
close MAIL;

# Return a new HTML page to the browser.  This page includes both the 
# download command and a thank-you page with links back to the software 
# package's home page.

print "Content-type:text/html\n\n";
print <<ENDHTML;

<html><head><title>Download Page</title></head>

<! The following line actually invokes the download >

<META HTTP-EQUIV="Refresh" CONTENT="0;URL=$downloadfile">
</HEAD>

<! The "thank you" page, using SNL template. >

<body text="#000000" background="http://www.sandia.gov/images/bkgrnd.gif">
<! KDD Turned off alternative link colors in template; the ><! following line was part of the above body command. ><! link="#003366" vlink="#cc0033" alink="#000000"><a NAME="TOP"></a><!---TOP BANNER AREA STARTS HERE--->
<table BORDER=0 valign="top" >
<tr VALIGN=TOP>
<td VALIGN=TOP WIDTH="140">
<table BORDER=0 WIDTH="130" valign="top" >
<tr VALIGN=TOP>
<td VALIGN=TOP WIDTH="128"><!--SANDIA LOGO AT TOP LEFT--><a href="http://www.sandia.gov/Main.html"><img SRC="http://www.sandia.gov/images/snlstkdc.gif" ALT="[Sandia National Laboratories]" BORDER=0 valign="top" height=49 width=126></a>
<p><img ISMAP SRC="http://www.sandia.gov/images/labelNEW.gif" ALT="[navigation panel]" HSPACE=2 BORDER=0 usemap="#shortMap" height=119 width=111></td>

<td><img SRC="http://www.sandia.gov/images/1pixel.gif" BORDER=0 height=1 width=10></td>
</tr>
</table>

<table BORDER=0 WIDTH="140" valign="top" >
<tr ALIGN=LEFT VALIGN=TOP>
<td VALIGN=TOP WIDTH="114"><!----------- 1st little turquoise bevel button ------------>
<table BORDER=0 CELLSPACING=0 CELLPADDING=0 WIDTH="114" BGCOLOR="#00CCFF" >
<tr VALIGN=TOP BGCOLOR="#99FFFF">
<td ALIGN=LEFT><img SRC="http://www.sandia.gov/images/1pixel.gif" height=1 width=1></td>

<td ALIGN=RIGHT><img SRC="http://www.sandia.gov/images/1pixel.gif" height=1 width=1></td>
</tr>

<tr ALIGN=CENTER VALIGN=CENTER>
<td COLSPAN="2"><b><font face="Verdana, Arial, Helvetica"><a href="$pkghome">
Return to $pkg pages</a></font></b></td>
</tr>

<tr VALIGN=BOTTOM BGCOLOR="#006699">
<td ALIGN=LEFT><img SRC="http://www.sandia.gov/images/1pixel.gif" height=1 width=1></td>

<td ALIGN=RIGHT><img SRC="http://www.sandia.gov/images/1pixel.gif" height=1 width=1></td>
</tr>
</table>
</td>

<td VALIGN=TOP WIDTH="20"></td>
</tr>

<tr VALIGN=TOP>
<td COLSPAN="2"></td>
</tr>
</table>
</td>

<td VALIGN=TOP><!--MAIN CONTENT AREA STARTS HERE--><!----------------THIS IS A CHANGE AREA----------------><!------HEADER TEXT SHOULD BE REPLACE THIS TEXT------><b><font face="Verdana, Arial, Helvetica"><font size=+2>Download Page&nbsp;</font></font></b>
<p><!---------------END OF THIS CHANGE AREA---------------><!----------------THIS IS A CHANGE AREA----------------><!--MAIN CONTENT SHOULD BE PLACED IN THE AREA BELOW-->
<hr width="100%">
<b>
Thank you for downloading 
<a href="$pkghome">$pkg</a>.
</b>
<p>

<b>To report downloading problems, contact:</b>
<ul>
<li>
<a href="mailto: kddevin\@cs.sandia.gov">Karen Devine</a></li>
</ul>
<!---------MAIN CONTENT AREA ENDS HERE---------><!-- CHANGE CONTACT + E-MAIL, NOTE "SUBJECT" IN E-MAIL CODE --></td>
</tr>
</table>

<table BORDER=0 WIDTH="100%" >
<tr ALIGN=CENTER>
<td VALIGN=TOP WIDTH="140">
<table BORDER=0 WIDTH="140" >
<tr>
<td ALIGN=CENTER VALIGN=TOP WIDTH="120"></td>

<td WIDTH="20"></td>
</tr>
</table>
</td>

<td ALIGN=CENTER VALIGN=TOP WIDTH="100%"></td>
</tr>
</table>
<!--Image maps below--><map name="shortMap"><area shape="rect" coords="2,2,108,14"href="http://www.sandia.gov/About.htm"></area><area shape="rect" coords="2,19,108,31"href="http://www.sandia.gov/Solution.htm"></area><area shape="rect" coords="2,36,108,48"href="http://www.sandia.gov/Working.htm"></area><area shape="rect" coords="2,53,108,65"href="http://www.sandia.gov/Contacting.htm"></area><area shape="rect" coords="2,70,108,82"href="http://www.sandia.gov/News.htm"></area><area shape="rect" coords="2,87,108,99"href="http://www.sandia.gov/search.html"></area><area shape="rect" coords="2,104,108,116"href="http://www.sandia.gov/Main.html"></area></map><!----------------THIS IS A CHANGE AREA----------------><!----NAME AND DATE OF LAST REVISION SHOULD BE HERE----><!---------------END OF THIS CHANGE AREA--------------->
</body>
</html>

ENDHTML

# end of Download page.


