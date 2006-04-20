<?php 

$counter = 1;

################################################################################

function print_header() {

  echo '<p class="heading">WebSolver</p>';
  
}

################################################################################

function print_problem_and_result($ProblemIDs, $ResultIDs, $flag)
{
  echo '<p class="subheading">Selected Problem IDs:</p>';
  echo '<ol>';
  $done = 0;
  foreach (explode(':', $ProblemIDs) as $i)
  {
    if ($i == "") continue;

    $j = explode('@', $i);
    if ($j[0] == "") continue;
    echo '<li>' . $j[0];
    $done = 1;
  }
  echo '</ol>';
  if ($done == 0)
    echo "<p>No problems are currently selected.";
  else
  {
    if ($flag == 1)
    { 
      echo '<form action="#" enctype="multipart/form-data" method="post" name="inputForm">';
      echo '<input type=hidden name=ProblemIDs value="">';
      echo '<input type=hidden name=ResultIDs value="' . $ResultIDs . '">';
      echo '<input type = submit class=submitPrimary value = "reset ProblemIDs" ></form>';
    }
  }

  echo '<p class="subheading">Recordered Results:</p>';
  echo '<ol>';
  $done = 0;
  foreach (explode(':', $ResultIDs) as $i)
  {
    if ($i == "") continue;

    $j = explode('@', $i);
    echo "<li>$j[1] &nbsp;&nbsp;-->&nbsp;&nbsp; phi = $j[0]";
    $done = 1;
  }
  echo '</ol>';
  if ($done == 0)
    echo "<p>No results are currently recordered.";
  else
  {
    if ($flag == 1)
    { 
      echo '<form action="#" enctype="multipart/form-data" method="post" name="inputForm">';
      echo '<input type=hidden name=ProblemIDs value="' . $ProblemIDs . '">';
      echo '<input type=hidden name=ResultIDs value="">';
      echo '<input type = submit class=submitPrimary value = "reset results" ></form>';
    }
  }
}

################################################################################

function fixed_parameter($name, $type, $value) 
{
  global $counter;
  echo '<tr>';
  $name2 = "name_" . $counter;
  echo "<td>" . $name . "</td>";
  echo "<input type=hidden name=$name2 value=\"$type:$name\"/></td>";
  $value2 = "value_" . $counter;
  echo "<td><input type=text name=$value2 value=\"$value\" size=50/></td>";
  echo '</tr>';
  $counter = $counter + 1;
}

################################################################################

function fixed_parameter2($name, $type, $value) 
{
  global $counter;
  $name2 = "name_" . $counter;
  echo "<input type=hidden name=$name2 value=\"$type:$name\"/></td>";
  $value2 = "value_" . $counter;
  echo "<td><input type=text name=$value2 value=\"$value\" size=10/></td>";
  $counter = $counter + 1;
}

################################################################################

function custom_parameter($name, $type, $value) 
{
  global $counter;
  echo '<tr>';
  $name2 = "name_" . $counter;
  echo "<td><input type=text name=$name2 value=\"$type:$name\" size=30/></td>";
  $value2 = "value_" . $counter;
  echo "<td><input type=text name=$value2 value=\"$value\" size=30/></td>";
  echo '</tr>';
  $counter = $counter + 1;
}

################################################################################

function process()
{ 
  $timestamp = date("y-m-d_H.i.s", time());

  global $ProblemIDs;

  $counter = $_POST['counter'];

  $configString  = "";
  $configString .= "ProblemIDs         = ".$ProblemIDs ."\n";
  $configString .= "i:iters            = ".$_POST['iters'] ."\n";
  $configString .= "d:tol              = ".$_POST['tol'] ."\n";
  $configString .= "s:az_solver        = ".$_POST['az_solver'] ."\n";
  $configString .= "b:perform_analysis = ".$_POST['perform_analysis'] ."\n";
  $configString .= "b:perform_cheby    = ".$_POST['perform_cheby'] ."\n";
  $configString .= "b:perform_jacobi   = ".$_POST['perform_jacobi'] ."\n";
  $configString .= "b:perform_gs       = ".$_POST['perform_gs'] ."\n";
  $configString .= "b:perform_sgs      = ".$_POST['perform_sgs'] ."\n";
  $configString .= "b:perform_ic       = ".$_POST['perform_ic'] ."\n";
  $configString .= "b:perform_ict      = ".$_POST['perform_ict'] ."\n";
  $configString .= "b:perform_ilu      = ".$_POST['perform_ilu'] ."\n";
  $configString .= "b:perform_ilut     = ".$_POST['perform_ilut'] ."\n";
  $configString .= "b:perform_ml       = ".$_POST['perform_ml'] ."\n";

  for ($i = 1; $i < $counter; $i++)
  {
    $name  = $_POST["name_"  . $i];
    $type  = $_POST["type_"  . $i];
    $value = $_POST["value_" . $i];
    $configString .= "$type$name = $value\n";
  }

  $configFile = fopen("/tmp/configs/$timestamp.txt", 'w')
    or die("can't open /tmp/configs/$timestamp.txt: $php_errormsg");
  if (-1 == fwrite($configFile, $configString)) { 
    die("can't write to /tmp/configs/$timestamp.txt: $php_errormsg"); }
  fclose($configFile) 
    or die("can't close /tmp/configs/$timestamp.txt: $php_errormsg");
  chmod("/tmp/configs/$timestamp.txt", 0664);

  chdir("solve/");

  $command = "PYTHONPATH=/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/:\$PYTHONPATH ";
  $command .= "python ./step_process.py /tmp/configs/$timestamp.txt 2>&1";
  passthru($command);
}


################################################################################

function step_header($thisID)
{ 
  global $ProblemIDs;
  global $ResultIDs;
?>

  <table>
  <tr><td>
  <form action="step1.html" enctype="multipart/form-data" method="post" name="inputForm">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 1: Data">
  </form> &nbsp; &nbsp;

  </td><td>

  <form action="#" enctype="multipart/form-data" method="post" name="step1">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="hidden" name="perform_analysis" value="True">
  <input type="submit" class=submitSecondary name="submit" value="Step 1C: Check Data">
  </form> &nbsp; &nbsp;

  </td><td>

  <form action="step2.html" enctype="multipart/form-data" method="post" name="step2">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 2: Parameters">
  </form> &nbsp; &nbsp;

  </td><td>

  <form action="#" enctype="multipart/form-data" method="post" name="step3">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 3: Compute">
  </form> &nbsp; &nbsp;

  </td><td>

  <form action="step4.html" enctype="multipart/form-data" method="post" name="step4">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 4: Results">
  </form> &nbsp; &nbsp;
  </td></tr></table>
<?
}
?>
