{a = csv_array("rect.csv")}
{_therow = rows(a)}
{_col = cols(a)}
{_r = 0}
{loop(_therow)}
{_c = 0}
{loop(_col)}
a[{_r},{_c}] = {a[_r,_c++]}
{endloop}
{_r++}
{endloop}
{print_array(a)}

{b = make_array(3,2)}
{NOECHO}
{_row = rows(b)}
{_col = cols(b)}
{_r = 0}
{loop(_row)}
{_c = 0}
{loop(_col)}
{b[_r,_c] = 10*(_r+1) + _c++ + 1}
{endloop}
{_r++}
{endloop}
{ECHO}
b = {print_array(b)}

{c = make_array(2,3)}
{NOECHO}
{_row = rows(c)}
{_col = cols(c)}
{_r = 0}
{loop(_row)}
{_c = 0}
{loop(_col)}
{c[_r,_c] = 10*(_r+1) + _c++ + 1}
{endloop}
{_r++}
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
{DUMP()}
