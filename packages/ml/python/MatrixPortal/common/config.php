<? 

# --------------------------------------------------- #
# File Locations                                      #
#                                                     #
# Specify here the location of the matrix to be       #
# uploaded and the script files.                      #
#                                                     #
# `MatrixDirectory' specifies where to store the      #
# uploaded matrices.                                  #
#                                                     #
# `ImageDirectory' specifies where to store the       #
# images created on-the-fly.                          #
#                                                     #
# `TempDirectory' is the directory used to store any  #
# temporary data.                                     #
#                                                     #
# `PythonDirectory' is the directory containing the   #
# step_process.py; typically, this is the ml/python/  #
# MatrixPortal directory.                             #
# --------------------------------------------------- #

$MatrixDirectory = "/var/www/html/MatrixPortal/HBMatrices/";
$ImageDirectory = "/var/www/html/tmp/";
$TempDirectory = "/tmp/";
$PythonDirectory = "/home/msala/Trilinos/packages/ml/python/MatrixPortal/";

# --------------------------------------------------- #
# Web Page Styles and Appearance                      #
#                                                     #
# Functions top_left() and top_right() can be used to #
# personalize the appearance of the web page.         #
# --------------------------------------------------- #

?>

<?  function top_left() { ?>
          <a href="<?php echo $dir ?>index.html"
          onmouseover="document.getElementById('trilinosTitle').src=trilinos_title_over.src;"
          onmouseout="document.getElementById('trilinosTitle').src=trilinos_title_normal.src;">
            <img class="trilinosTitle" id="trilinosTitle" src="<?php echo $dir ?>common/trilinos_title_normal.png"
            alt="The Trilinos Project" /></a>
<?  } ?>

<? function top_right() { ?>
        <a href="http://www.sandia.gov/" target="_blank"
          onmouseover="document.getElementById('sandia').src=sandia_over.src;"
          onmouseout="document.getElementById('sandia').src=sandia_normal.src;">
            <img class="sandia" id="sandia" src="<?php echo $dir ?>common/sandia_normal.png"
            alt="Sandia National Laboratories" /></a>
<?  } ?>
