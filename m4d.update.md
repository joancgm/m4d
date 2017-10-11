__update October 10, 2017__
__corrections to m4d.2015.2 --> version m4d.2015.2c__

__Corrections - code__
cvdcsubs.c (c_cvdcinit), change dimension of i4dp from 3 to 4   
c_set_cpda.c  change dimension of idd and idpc from 3 to 4      
c_varinit.c correct type 1 error when variable is not new 
c_varinit.c fix so input file closes for types 2 and 3 

__Corrections - documentation__
m4d.commands.doc (and .pdf) - equation for cvdci in cvdcreset
m4d.commands.short.doc (and .pdf)  lineplot input  imageave before lineave

__Modification__
c_varinit.c   change readname to readfilename so can use string alias for types 2 and 3
m4d.c  label as version 2015.2c in printout
m4d.commands.doc (and .pdf) corresponding corrections and changes


