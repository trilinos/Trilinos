 GENERATE A PLOT.
** CHECKOUT.   
 LEGEND OFF.
 PAGE-BORDER OFF.
 X TICKS 5.
 Y TICKS 5.
 EVERY CURVE SYMBOL COUNT 0.
 CURVE 1 TEXTURE 1 COLOR RED THICKNESS 5.
 CURVE 2 TEXTURE 1 COLOR GREEN thickness 5.
 CURVE 3 symbol count 1, COLOR yellow.
 CURVE 4 TEXTURE dotted COLOR GREEN.
 CURVE 5 TEXTURE dotted COLOR GREEN.
 EVERY MESSAGE UNITS COORDINATE, BLANKING OFF.
 message 1 text "" x   -.03   y  -0.1.
 MESSAGE 1 pointer cc -0.2297 0.9732 ARROWHEAD 0401 UNITS COORDINATE.
 message 2 text "S2" x -0.25 y 0.35 connect tr.
 message 3 text "S1" x  0.25 y 0.15 connect bl.
 message 4 text "M"  x  0.0 y -0.1 connect tc.
 message 5 text "Normal Vector" x -0.15 y 0.9 connect lc.
 message 6 text "Master Surface" x 0.1  y -0.15 connect tl.
 message 7 text "Slave Surface" x 0.1   y 0.25 connect bl.
 AXIS CROSS OFF.
 X AXIS MINIMUM  -1.0000E+00, MAXIMUM  1.
 Y AXIS MINIMUM  -1.0000E+00, MAXIMUM  1.
 grid color cyan.
 INPUT DATA.
"CURVE 1"
-0.5  -0.35
 0.0  -0.1
 0.5  -0.1
"2"
-0.5   0.5
 0.5   0.0
"3"
 -0.2 0.35
  0.2 0.15
"4" 
  0.0  -0.1
  0.2   0.15
END OF DATA.
