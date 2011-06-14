    static char *qainfo[] = 
      {
       "Aprepro ",
       "08/07/90",
       "X1.01.04"
      };
/* x1.01.01: 07/02/90  16:03:37  GREG SJAARDEMA2                            */
/*         : Initial installation of automatic QA updating in Aprepro       */
/* X1.01.02: 07/10/90  15:14:42  GREG SJAARDEMA1                            */
/*         : Added ECHO/NOECHO commands to control output to output file.   */
/*         : ECHO causes all lines to be echoed to output file (default)    */
/*         : NOECHO turns off echoing, lines are still processed.           */
/* X1.01.03: 07/11/90  08:39:32  GREG SJAARDEMA2                            */
/*         : Modified ECHO/NOECHO slightly to use True/False instead of     */
/*         :  1/0, Turn ECHO back on at end of included files.              */
/* X1.01.04: 08/07/90  10:25:47  GREG SJAARDEMA3                            */
/*         : Added VERBATIM(ON|OFF) command. Controls processing of input   */
/*         : lines.  VERBATIM(ON) will echo all lines to output without     */
/*         : processing.  VERBATIM(OFF) returns to normal processing.       */
