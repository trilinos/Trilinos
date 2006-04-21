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
            <a href="index.html" target="_blank" 
            <img src=common/matrix_portal_logo.png height=60 border=0 alt="The Matrix Portal" /></a>
<?  } ?>

<? function top_right() { ?>
        <a href="http://www.ethz.ch/" target="_blank" 
            <img src=common/eth_logo.gif border=0 alt="ETHZ" /></a>
<?  } ?>
