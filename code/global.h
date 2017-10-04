#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

int printout(const char *type, const char *format, ...);  /* in m4d.c */
void exitm4d(int i);   /* in m4d.c */
  
#include "arrays.h"  /* array handling routines */
#include "input.h"   /* read routines */
#include "interp.h"  /* interpolation routines */
#include "geom.h"    /* geom8 and geomc routines */
#include "vector.h"  /* extra math definitions     */
#include "coefsubs.h" /* non-command coef subs */
#include "iexpand.h"
#include "fixdudx.h"

