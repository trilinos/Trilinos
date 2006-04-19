<?php 

$counter = 1;

################################################################################

function print_header() {

  echo '<p class="heading">WebSolver</p>';
  
}

################################################################################

function print_problem($ProblemIDs, $flag) {

  echo '<p class="subheading">Current Settings:</p>';
  echo 'Selected problem IDs:';
  echo '<ol>';
  foreach (explode(':', $ProblemIDs) as $i)
  {
    if ($i == "") continue;
    echo '<li>' . $i;
  }
  echo '</ol>';
  if ($flag == 1)
  { 
    echo '<form action="#" enctype="multipart/form-data" method="post" name="inputForm">';
    echo '<input type=hidden type=ProblemIDs value="">';
    echo '<input type = submit value = "reset ProblemIDs" ></form>';
  }
  
}

################################################################################

function fixed_parameter($name, $type, $value) {

  global $counter;
  echo '<tr>';
  $name2 = "__PyTrilinos__name" . $counter;
  echo "<td>" . $name . "</td>";
  echo "<input type=hidden name=$name2 value=\"$type:$name\"/></td>";
  $value2 = "__PyTrilinos__value" . $counter;
  echo "<td><input type=text name=$value2 value=\"$value\" size=50/></td>";
  echo '</tr>';
  $counter = $counter + 1;

}

################################################################################

function fixed_parameter2($name, $type, $value) {

  global $counter;
  $name2 = "__PyTrilinos__name" . $counter;
  echo "<input type=hidden name=$name2 value=\"$type:$name\"/></td>";
  $value2 = "__PyTrilinos__value" . $counter;
  echo "<td><input type=text name=$value2 value=\"$value\" size=10/></td>";
  $counter = $counter + 1;

}

################################################################################

function custom_parameter($name, $type, $value) {

  global $counter;
  echo '<tr>';
  $name2 = "__PyTrilinos__name" . $counter;
  echo "<td><input type=text name=$name2 value=\"$type:$name\" size=30/></td>";
  $value2 = "__PyTrilinos__value" . $counter;
  echo "<td><input type=text name=$value2 value=\"$value\" size=30/></td>";
  echo '</tr>';
  $counter = $counter + 1;

}

################################################################################

function process($analyze, $direct, $iterative) { 
  
  chdir("/people_old/trilinos_www/htdocs/");

  $timestamp = date("y-m-d_H.i.s", time());

  global $ProblemIDs;

  $counter = $_POST['counter'];
  $configString = "";
  $configString .= "COUNTER = " . $counter . "\n";
  $configString .= "SUBMIT_TIME = ".$timestamp."\n";
  $configString .= "EMAIL = \n";
  $configString .= "\n";
  $configString .= "PROBLEM_ID = ".$ProblemIDs ."\n";
  $configString .= "DIRECT_SOLVER = ".$_POST['direct_solver'] ."\n";
  $configString .= "ITERATIVE_SOLVER = ".$_POST['iterative_solver'] ."\n";
  $configString .= "ITERS = ".$_POST['iters'] ."\n";
  $configString .= "TOL = ".$_POST['tol'] ."\n";

  for ($i = 1; $i < $counter; $i++)
  {
    $name = "__PyTrilinos__name" . $i;
    $configString .= "$name = ".$_POST[$name] ."\n";
    $type = "__PyTrilinos__type" . $i;
    $configString .= "$type = ".$_POST[$type] ."\n";
    $value = "__PyTrilinos__value" . $i;
    $configString .= "$value = ".$_POST[$value] ."\n";
  }

  $configFile = fopen("/tmp/configs/$timestamp.txt", 'w')
    or die("can't open /tmp/configs/$timestamp.txt: $php_errormsg");
  if (-1 == fwrite($configFile, $configString)) { 
    die("can't write to /tmp/configs/$timestamp.txt: $php_errormsg"); }
  fclose($configFile) 
    or die("can't close /tmp/configs/$timestamp.txt: $php_errormsg");
  chmod("../configs/$timestamp.txt", 0664);

  chdir("solve/");

  if ($analyze == 1)
  {
    echo "<div class=\"outputBox\"><pre>";
    $command = "PYTHONPATH=/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/:\$PYTHONPATH ";
    $command .= "python ./step_process.py /tmp/configs/$timestamp.txt analyze 2>&1";
    passthru($command);
    echo "&nbsp;<pre></div>";
  }

  if ($direct == 1)
  {
    echo "<div class=\"outputBox\"><pre>";
    $command = "PYTHONPATH=/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/:\$PYTHONPATH ";
    $command .= "python step_process.py /tmp/configs/$timestamp.txt direct 2>&1";
    passthru($command);
    echo "&nbsp;<pre></div>";
  }

  if ($iterative == 1)
  { 
    echo "<div class=\"outputBox\"><pre>";
    $command = "PYTHONPATH=/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/:\$PYTHONPATH ";
    $command .= "python step_process.py /tmp/configs/$timestamp.txt iterative 2>&1";
    passthru($command);
    echo "&nbsp;<pre></div>";
  }
  
  chdir("/people_old/trilinos_www/htdocs/");

}

################################################################################

function my_footer() {

  global $ProblemIDs;

  echo "<br>";
  echo "<P>";
  echo "<table>";
  echo '<tr><td><form action="step_0.html" enctype="multipart/form-data" method="post" name="inputForm">';
  echo '<input type=submit value="HOME" ></form></td></tr>';
  
  echo '<tr><td><form action="step_1.html" enctype="multipart/form-data" method="post" name="inputForm">';
  echo '<input type=hidden name=ProblemIDs value="' . $ProblemIDs . '">';
  echo '<input type=submit value="Step 1" ></form></td></tr>';

}

################################################################################

function write($what) {

  global $ProblemID;

  $filename = '/people_old/trilinos_www/matrices/comments.txt';
  $fp = fopen($filename, "a");
  fputs($fp, '<p><b>' . $ProblemID . '</b>');
  fputs($fp, '<p>' . $what);
  fclose($fp);

}

################################################################################

function comment() {

  global $ProblemID;
  
  ?>
  <form action="#" enctype="multipart/form-data" method="post" name="writeComment">

  <p>Comments to be added to step 4:
  <input type="hidden" name="stage" value="writeComment">
  <input type="hidden" name="problemID" value="<? global $ProblemID; echo $ProblemID; ?>" />
  <p><textarea name=comment rows=4 cols=60></textarea>
  <input type=submit value="add to comments">
  </form>
  <?
  
} 

?>
