#include <stdlib.h>
int sierra_xmldb_original_length = 1 ;
static unsigned char sierra_xmldb[] = {(unsigned char)0xff};
size_t sierra_xmldb_buffer_size() {return sizeof(sierra_xmldb)/sizeof(sierra_xmldb[0]);}
unsigned char*  sierra_xmldb_buffer() { return &sierra_xmldb[0]; }
