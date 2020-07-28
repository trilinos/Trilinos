/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* devid.h - these functions are used to map SVDI device code numbers
 *           to device code character strings and vice versa
 * Debbie Campbell
 */
struct device_ids
{
  char *devid_char;
  float devid_num;
};

/*  to add new devices to this table:  increment MAX_DEVID value
                                       add to table per example below    */
#define MAX_DEVID 74
static struct device_ids device_values[MAX_DEVID] = {
    {"tk4", 1.},   {"tk6", 1.1},  {"tek", 1.2},  {"tp2", 1.3},  {"tp8", 1.4},  {"t14", 1.5},
    {"t13", 1.6},  {"t05", 1.7},  {"t07", 1.8},  {"t15", 1.9},  {"16c", 2.1},  {"16b", 2.2},
    {"35c", 2.3},  {"3mc", 2.31}, {"810", 2.32}, {"35b", 2.4},  {"3mb", 2.41}, {"35a", 2.5},
    {"24l", 2.6},  {"48l", 2.7},  {"csq", 2.8},  {"bsq", 2.9},  {"r94", 3.},   {"t27", 4.},
    {"t25", 4.1},  {"alp", 5.},   {"hc1", 6.},   {"lxy", 7.},   {"tst", 8.},   {"v25", 9.},
    {"v40", 9.1},  {"v41", 9.2},  {"v34", 9.21}, {"aed", 10.},  {"ae7", 10.1}, {"ae1", 10.2},
    {"met", 11.},  {"hpp", 12.},  {"h75", 12.1}, {"h72", 12.2}, {"h74", 12.3}, {"h50", 12.4},
    {"ret", 14.},  {"ap5", 15.},  {"jp7", 16.},  {"jp1", 16.1}, {"ger", 17.},  {"xyn", 18.},
    {"ps3", 20.},  {"qms", 21.},  {"c51", 22.},  {"c10", 22.1}, {"r25", 23.},  {"r38", 23.1},
    {"f3c", 23.2}, {"f3m", 23.3}, {"f16", 23.4}, {"fsq", 23.5}, {"f8t", 23.6}, {"qlf", 24.},
    {"q35", 24.1}, {"t45", 25.},  {"ls5", 26.},  {"i10", 27.},  {"i30", 27.1}, {"btk", 28.},
    {"ln3", 29.},  {"a60", 30.},  {"sun", 31.},  {"uis", 31.1}, {"ult", 32.},  {"tal", 33.},
    {"nws", 34.},  {"x10", 35.}};

#if !defined(NO_GET_DEVID_CHAR)
/*                                                                          *
 *  function get_devid_char will return a pointer to the character string   *
 *  associated with a floating point number for cgid inquiries.             *
 *               INPUT: a pointer to the 3 character string                 *
 *               RETURN: floating point device number                       *
 *                                                                          */
static char *get_devid_char(float number)
{
  int i;

  for (i = 0; i < MAX_DEVID; i++) {
    if (device_values[i].devid_num == number) {
      return (device_values[i].devid_char);
    }
  }

  /*  return a null pointer if there is no floating point number match     */
  return (0);
}
#endif
/* end devid.h */
