{a = csv_array("rect.csv")}
{_therow = rows(a)}
{_col = cols(a)}
{loop(_therow, _r)}
{loop(_col, _c)}
a[{_r},{_c}] = {a[_r,_c]}
{endloop}
{endloop}
{print_array(a)}

{b = make_array(3,2)}
{NOECHO}
{_therow = rows(b)}
{_col = cols(b)}
{loop(_therow, _r)}
{loop(_col, _c)}
{b[_r,_c] = 10*(_r+1) + _c + 1}
{endloop}
{endloop}
{ECHO}
b = {print_array(b)}

{c = make_array(2,3)}
{NOECHO}
{_therow = rows(c)}
{_col = cols(c)}
{loop(_therow, _r)}
{loop(_col, _c)}
{c[_r,_c] = 10*(_r+1) + _c + 1}
{endloop}
{endloop}
{ECHO}
c = {print_array(c)}

{d = b * c}
d = b * c = {print_array(d)}

{e = c * b}

e = c * b = {print_array(e)}
{et = transpose(e)}

e' = {print_array(transpose(e))}

sum = e + e = {print_array(e+e)}

-e = {print_array(-e)}

e+e-2*e = {print_array(e+e-2*e)}

scale = {print_array(e*1/e[0,0])}

{asub = csv_array("rect.csv", 10)}
{asub_rows = rows(asub)}
{a_rows = rows(a)}

{std = array_from_string('2.3 5.7 11.13 17.19 23.29 31.37 41.43 47.53 59.61 67.71',' ')}
{print_array(transpose(std))}

{st = array_from_string('2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71',' ')}
{print_array(transpose(st))}

{DUMP()}
