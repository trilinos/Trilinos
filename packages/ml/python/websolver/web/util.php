<!-- ####################################################################### -->

<?php function connectToDatabase() { 

  // open the database connection
  if (!($connection = @ mysql_connect("localhost","websolver"))) {
    die ("Error: unable to connect to database.");
  }
  
  // select the trilinos database
  if (!(@ mysql_select_db("trilinos", $connection))) {
    printError();
  }
       
  return $connection;
      
} ?>

<!-- ####################################################################### -->

<?php function printError() { ?>

  <span style="font-size: 120%; font-weight: bold; color: #f00;">
    MySQL Error #<?php echo mysql_errno(); ?> <br /> <br />
    <?php echo mysql_error(); ?>
  </span>
  
<?php } ?>