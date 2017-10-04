/* contains c_comment */
#include "global.h"

/* read and write a one line comment */
void c_comment(FILE *fpin, FILE *fprint)
{ 
  char *comment;
  comment=readline(fpin);
  printout("print","%s\n",comment);
}
