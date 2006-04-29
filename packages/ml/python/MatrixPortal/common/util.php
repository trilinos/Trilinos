<?php 

$counter = 1;

################################################################################

function print_header() {

  echo '<p class="heading">WebSolver</p>';
  
}

################################################################################

function print_problem_and_result($ProblemIDs, $ResultIDs, $flag)
{
?>
  <p class'subheading'><a onclick="toggle('IDs');" onkeypress="toggle('IDs');"><b>Labels and Results</b></p></a>
<?
  echo '<div class=open id=IDs>';
  echo '<form action="#" enctype="multipart/form-data" method="post">';
  echo '<input type=hidden name=ProblemIDs value="' .  $ProblemIDs . '">';
  echo '<input type=hidden name=ResultIDs value="' .  $ResultIDs . '">';
  echo '<table border=0><cols=3>';
  echo '<tr><td>Problem IDs</td>';
  echo '<td>&nbsp;&nbsp;&nbsp;</td>';
  echo '<td>Recordered Results</td></tr>';
  echo '<td valign=top><ol>';
  $done = 0;
  $done2 = 0;
  $count = 1;
  foreach (explode(':', $ProblemIDs) as $i)
  {
    if ($i == "") continue;

    echo '<li><input type=checkbox name="dp_' . $count++ . '" value=Yes>&nbsp;';

    $j = explode('@', $i);
    if ($j[0] == "") continue;
    echo $j[0];
    $done = 1;
  }

  echo '</ol>';
  $count = 1;
  if ($done == 0)
	  echo "No problems are currently selected.";

	  echo '</td><td></td><td valign=top>';
	  echo '<ol>';
	  foreach (explode(':', $ResultIDs) as $i)
	  {
	    if ($i == "") continue;

            echo '<li><input type=checkbox name="dr_' . $count++ . '" value=Yes>&nbsp;';

	    $j = explode('@', $i);
	    echo "$j[1] &nbsp;&nbsp;<font color=red>phi = $j[0]</font>";
	    $done2 = 1;
	  }
	  echo '</ol>';
	  if ($done2 == 0)
	    echo "No results are currently recordered.";

  echo '</td></tr></table>';
  if ($done != 0 || $done2 != 0)
  {
    echo '<table border=0><tr valign=top><td>';
    echo '<input type=submit class=submitPrimary value="delete selected"></form>';
    echo '</td><td>';
    echo '<form action="#" enctype="multipart/form-data" method="post" name="inputForm">';
    echo '<input type=hidden name=ProblemIDs value="">';
    echo '<input type=hidden name=ResultIDs value="' . $ResultIDs . '">';
    echo '<input type = submit class=submitPrimary value = "reset all ProblemIDs" ></form>';
    echo '</td><td>';
    echo '<form action="#" enctype="multipart/form-data" method="post" name="inputForm">';
    echo '<input type=hidden name=ProblemIDs value="' . $ProblemIDs . '">';
    echo '<input type=hidden name=ResultIDs value=>';
    echo '<input type = submit class=submitPrimary value = "reset all ResultIDs" ></form>';
    echo '</td></tr></table>';
  }
  echo '</div>';
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

function begin_select_parameter($info, $name, $type)
{
  global $counter;
  echo '<tr>';
  $name2 = "name_" . $counter;
  echo "<td>" . $info . "</td>";
  echo "<input type=hidden name=$name2 value=\"$type:$name\"/></td>";
  $value2 = "value_" . $counter;
  echo "<td><select name=\"" . $value2 . "\">";
}

################################################################################

function add_select_parameter($value, $desc) 
{
  global $counter;
  $value2 = "value_" . $counter;
  echo '<option value="' . $value . '">' . $desc;
  $counter = $counter + 1;
}

################################################################################

function end_select_parameter()
{
  echo '</select></td></tr>';
}

################################################################################

function begin_select_parameter_notable($name, $type)
{
  global $counter;
  $name2 = "name_" . $counter;
  echo "<input type=hidden name=$name2 value=\"$type:$name\"/>";
  $value2 = "value_" . $counter;
  echo "<select name=\"" . $value2 . "\">";
}

################################################################################

function add_select_parameter_notable($value, $desc) 
{
  global $counter;
  $value2 = "value_" . $counter;
  echo '<option value="' . $value . '">' . $desc;
  $counter = $counter + 1;
}

################################################################################

function end_select_parameter_notable()
{
  echo '</select>';
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
  global $TempDirectory;
  global $PythonDirectory;

  $counter = $_POST['counter'];

  $configString  = "";
  $configString .= "ProblemIDs         := ".$ProblemIDs ."\n";
  $configString .= "i:iters            := ".$_POST['iters'] ."\n";
  $configString .= "d:tol              := ".$_POST['tol'] ."\n";
  $configString .= "s:az_solver        := ".$_POST['az_solver'] ."\n";
  $configString .= "i:az_kspace        := ".$_POST['az_kspace'] ."\n";
  $configString .= "s:az_output        := ".$_POST['az_output'] ."\n";
  $configString .= "s:rhs              := ".$_POST['rhs'] ."\n";
  $configString .= "s:starting_solution:= ".$_POST['starting_solution'] ."\n";
  $configString .= "s:solution         := ".$_POST['solution'] ."\n";
  $configString .= "b:perform_analysis := ".$_POST['perform_analysis'] ."\n";
  $configString .= "b:perform_cheby    := ".$_POST['perform_cheby'] ."\n";
  $configString .= "b:perform_jacobi   := ".$_POST['perform_jacobi'] ."\n";
  $configString .= "b:perform_gs       := ".$_POST['perform_gs'] ."\n";
  $configString .= "b:perform_sgs      := ".$_POST['perform_sgs'] ."\n";
  $configString .= "b:perform_ic       := ".$_POST['perform_ic'] ."\n";
  $configString .= "b:perform_ict      := ".$_POST['perform_ict'] ."\n";
  $configString .= "b:perform_ilu      := ".$_POST['perform_ilu'] ."\n";
  $configString .= "b:perform_ilut     := ".$_POST['perform_ilut'] ."\n";
  $configString .= "b:perform_ml       := ".$_POST['perform_ml'] ."\n";

  for ($i = 1; $i < $counter; $i++)
  {
    $name  = $_POST["name_"  . $i];
    $type  = $_POST["type_"  . $i];
    $value = $_POST["value_" . $i];
    if ("$type$name" != "")
      $configString .= "$type$name := $value\n";
  }

  $configFileName = "$TempDirectory/configs/$timestamp.txt";
  $configFile = fopen($configFileName, 'w')
    or die("can't open $configFileName: $php_errormsg");
  if (-1 == fwrite($configFile, $configString)) { 
    die("can't write to $configFileName: $php_errormsg"); }
  fclose($configFile) 
    or die("can't close $configFileName: $php_errormsg");
  chmod($configFileName, 0664);

  $command = "PYTHONPATH=/home/msala/Trilinos/LINUX_SERIAL/lib/python2.4/site-packages/:\$PYTHONPATH ";
  $command .= "python $PythonDirectory/step_process.py $configFileName 2>&1";
  passthru($command);
}


################################################################################

function step_header($thisID)
{ 
  global $ProblemIDs;
  global $ResultIDs;
?>

  <table border=0>
  <tr><td>
  <form action="step1.html" enctype="multipart/form-data" method="post" name="inputForm">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 1: Select Data">
  </form> &nbsp; &nbsp;

  </td><td>


  <form action="step2.html" enctype="multipart/form-data" method="post" name="step2">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 2: Select Parameters">
  </form> &nbsp; &nbsp;

  </td><td>

  <? if ($thisID != "3") { ?>
  <form action="step4.html" enctype="multipart/form-data" method="post" name="step4">
  <input type="hidden" name="ProblemIDs" value="<? global $ProblemIDs; echo $ProblemIDs; ?>">
  <input type="hidden" name="ResultIDs" value="<? global $ResultIDs; echo $ResultIDs; ?>">
  <input type="submit" class=submitSecondary name="submit" value="Step 4: Check Results">
  </form> &nbsp; &nbsp;
  <? } ?>
  </td></tr></table>
<?
}

################################################################################

function process_variables()
{ 
  global $ProblemIDs;
  global $ResultIDs;
  $OldProblemIDs = $_POST['ProblemIDs'];
  $OldResultIDs = $_POST['ResultIDs'];
  
  $ProblemIDs = "";
  $ResultIDs = "";

  $count = 1;
  foreach (explode(':', $OldProblemIDs) as $i)
  {
    if ($i == "") continue;
    if ($_POST['dp_' . $count++] != "Yes")
      $ProblemIDs .= ":" . $i;
  }

  $count = 1;
  foreach (explode(':', $OldResultIDs) as $i)
  {
    if ($i == "") continue;
    if ($_POST['dr_' . $count++] != "Yes")
      $ResultIDs .= ":" . $i;
  }
}

################################################################################

function check_data()
{ 
  global $ProblemIDs;
  global $ResultIDs;

  $OldProblemIDs = $ProblemIDs;
  $OldResultIDs = $ResultIDs;

  $ProblemIDs = "";

  foreach (explode(':', $OldProblemIDs) as $i)
  {
    if ($i == "") continue;
    $j = explode('@', $i);
    if (strstr($j[0], ":") ||
        strstr($j[0], "@") ||
        strstr($j[0], ";") ||
        strstr($j[0], "!") ||
        strstr($j[0], "+") ||
        strstr($j[0], "\\") ||
        strstr($j[0], "_") ||
        strstr($j[0], ":="))
    {
      echo "<p><font color=red><b>Sorry, ProblemIDs cannot contain \ : ; _ +  @ :=";
      echo ". The ProblemID was " . $j[0] . "</b></font></p>";
    }
    else
      $ProblemIDs .= $i . ":";
  }

  $ResultIDs = "";

  foreach (explode(':', $OldResultIDs) as $i)
  {
    if ($i == "") continue;
    $j = explode('@', $i);
    if (strstr($j[1], ":") ||
        strstr($j[1], "@") ||
        strstr($j[1], ";") ||
        strstr($j[1], "!") ||
        strstr($j[1], "+") ||
        strstr($j[1], "\\") ||
        strstr($j[1], "_") ||
        strstr($j[1], ":="))
    {
      echo "<p><font color=red><b>Sorry, ResultIDs cannot contain \ : ; _ +  @ :=";
      echo ". The resultID was " . $j[1] . "</b></font></p>";
    }
    else
      $ResultIDs .= $i . ":";
  }

}
?>
