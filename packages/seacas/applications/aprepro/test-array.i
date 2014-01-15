{_FORMAT = "%.3f"}
{a = csv_array("rect.csv")}
{row = rows(a)}
{col = cols(a)}
{_r = 0}
{loop(row)}
{_c = 0}
{loop(col)}
a[{_r},{_c}] = {a[_r,_c++]}
{endloop}
{_r++}
{endloop}
{print_array(a)}

{b = make_array(3,2)}
{row = rows(b)}
{col = cols(b)}
{_r = 0}
{loop(row)}
{_c = 0}
{loop(col)}
{b[_r,_c] = 10*(_r+1) + _c++ + 1}
{endloop}
{_r++}
{endloop}
b = {print_array(b)}

{c = make_array(2,3)}
{row = rows(c)}
{col = cols(c)}
{_r = 0}
{loop(row)}
{_c = 0}
{loop(col)}
{c[_r,_c] = 10*(_r+1) + _c++ + 1}
{endloop}
{_r++}
{endloop}
c = {print_array(c)}

{d = b * c}
d = {print_array(d)}

{e = c * b}

e = {print_array(e)}
{et = transpose(e)}

e' = {print_array(transpose(e))}

sum = {print_array(e+e)}
-e = {print_array(-e)}

scale = {print_array(e*1/e[0,0])}
{DUMP()}
