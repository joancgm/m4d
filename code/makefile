#  Choose either gcc or g++ compile by commenting in/out lines with #
#  These are incompatible so must do 'make clean' when changing.
#   October 2017, needs changes to compile under c++, the compiler objects
#   to comparing  a pointer (double *) to the integer 0.
#  If compiler doesn't understand CFLAGS, used reduced set of options.

PROGRAM       = ../a.m4d

#CC	      = g++
#LD	      = g++

CC	      = gcc
LD	      = gcc

CFLAGS	      = -Wall -O3 -Wuninitialized -ffast-math

#CFLAGS	      = -Wall -O3

DEST	      = .

EXTHDRS	      =

HDRS	      = arrays.h blocksubs.h commandlist.h cmatchfollow.h coefsubs.h \
                eqnsubs.h fixbij.h fixdudx.h geom.h global.h iexpand.h input.h\
                interp.h psleepcalc.h psolve1d.h tdmasubs.h vector.h p_plot.h

LDFLAGS	      = 

LIBS	      = -lm

MAKEFILE      = makefile

OBJS	      = m4d.o \
                commandlist.o \
                arrays.o \
                blocksubs.o \
		cmatchfollow.o \
		coefsubs.o \
		cvdcsubs.o \
		eqnsubs.o \
		exitsubs.o \
		fixbij.o \
		fixdudx.o \
		geom8subs.o \
		geomcsubs.o \
		iexpand.o \
		inletsubs.o \
		input.o \
		interp.o \
		omwallsubs.o \
		prinstress.o \
		psleepcalc.o \
		psolve1d.o \
		readmefpgrid.o \
		tdmasubs.o \
		c_algebra.o \
		c_areaflowint.o \
		c_arraydelete.o \
		c_arraydump.o \
		c_arrayhowto.o \
		c_arrayread.o \
		c_aveijk.o \
		c_bijrhsmarv.o \
		c_bijwallmarv.o \
		c_bijwallsplat.o \
		c_coefadd.o \
		c_coefconv.o \
		c_coefcplus.o \
		c_coefdt.o \
		c_coeffix.o \
		c_coefinit.o \
		c_coefrhs.o \
		c_coefvisc.o \
		c_coefzero.o \
		c_comment.o \
		c_constant.o \
		c_contcpcdu.o \
		c_contcpcfixed.o \
		c_contdu.o \
		c_contrhsp.o \
		c_contrhsu.o \
		c_copy.o \
		c_copypart.o \
		c_copy_mtop.o \
		c_copy_ptom.o \
		c_ddtall.o \
		c_edit.o \
		c_editabcd.o \
		c_eqpts2gpts.o \
		c_eqppts2ppts.o \
		c_eqnsolvebij.o \
		c_eqnsolvep.o \
		c_eqnsolves.o \
		c_eqnupdatem.o \
		c_function.o \
		c_gpts2eqpts.o \
		c_gradprop.o \
		c_gridcorner.o \
		c_gridfrommefp.o \
		c_gridmatch.o \
		c_gridtomefp.o \
		c_gridvarmod.o \
		c_if.o \
		c_interp_gtom.o \
		c_interp_gtop.o \
		c_interp_ptod.o \
		c_interp_ptog.o \
		c_lineoutput.o \
		c_lineoutputijk.o \
		c_mirror.o \
		c_momcam.o \
		c_momcamddt.o \
		c_momrhsr.o \
		c_omrhsmarv.o \
		c_omrhscoakley.o \
		c_ppreset.o \
		c_ppts2eqppts.o \
		c_prinstress.o \
		c_print.o \
		c_prints.o \
		c_printscoef.o \
		c_qturbrhs.o \
		c_rmsminmax.o \
		c_set_cpda.o \
		c_set_cpflop.o \
		c_set_cpsleep.o \
		c_set_gbij.o \
		c_set_pkdk.o \
		c_set_srate.o \
		c_set_volcont.o \
		c_set_volmom.o \
		c_set_wherefw.o \
		c_set_wherep.o \
		c_set_xyzdouble.o \
		c_valueat.o \
		c_varinit.o \
		c_varmatch.o \
		c_varupdate.o \
		c_vint2eqpts.o \
		c_visccoakley.o \
		c_viscmarv.o \
		c_viscmarvheat.o \
		c_viscqnoise.o \
		c_walldist.o \
		c_wallflux.o \
		c_wallnorm.o \
		ec_bijrhsmarvex.o \
		ec_qbijles.o \
		ec_rayleigh.o \
		ec_viscles.o \
		ec_w2wlineflat.o \
		pc_abcmask.o \
		pc_bar.o \
		pc_giflist.o \
		pc_gridinfo.o \
		pc_image.o \
		pc_iplaneint.o \
		pc_keyword.o \
		pc_lineplot.o \
		pc_picture.o \
		pc_xytocgrid.o \
		pc_xyzicut.o \
		p_fillcontour.o \
		p_glines.o \
		p_nextpt2d.o \
		p_outgif.o \
		p_plotsubs.o \
		p_symbol.o \
		p_vectors.o
           
PRINT	      = pr

#PROGRAM       = ../a.m4d

SHELL	      = /bin/sh

SRCS	      = 

SYSHDRS	      = 

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS)  -o $(PROGRAM) $(LIBS)
		echo "done"

DEPFILES=$(subst .o,.d,$(OBJS)) 

$(DEPFILES) : $(HDRS)

clean:
	rm $(OBJS) $(DEPFILES)

# The following creates the dependency, .d files

%.d : %.c
	$(CC) -MM $< | sed -e "s///g" > $@

include $(DEPFILES)



