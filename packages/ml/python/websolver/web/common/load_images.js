<!--// load_images.js

/* *****************************************************************************

    - Pre-loads rollover images into local cache
    - Usage (place in head section):
        <script src="load_images.js" language="JavaScript" type="text/javascript">
        </script>
    - Works the same in Explorer (tested only on 6.0) and Netscape (tested only
      on Netscape 7.0 and Mozilla 1.3.1), but does not work in Opera (6.01)

***************************************************************************** */

if (document.images) { // tests for image support
  
  trilinos_title_normal = new Image(310, 55);
  trilinos_title_normal.src = "<?php echo $dir ?>common/trilinos_title_normal.png";
  trilinos_title_over = new Image(310, 55);
  trilinos_title_over.src = "<?php echo $dir ?>common/trilinos_title_over.png";

  sandia_normal = new Image(180, 29);
  sandia_normal.src = "<?php echo $dir ?>common/sandia_normal.png";
  sandia_over = new Image(180, 29);
  sandia_over.src = "<?php echo $dir ?>common/sandia_over.png";
  
  trilinos_normal = new Image(183, 100);
  trilinos_normal.src = "<?php echo $dir ?>common/trilinos_normal.png";
  trilinos_over = new Image(183, 100);
  trilinos_over.src = "<?php echo $dir ?>common/trilinos_over.png";
  
} // if (document.images)

//-->
