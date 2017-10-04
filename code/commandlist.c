/* contains commandlist */
#include "global.h"
#include "commandlist.h"   /* list of commands all type void */
#include "p_plot.h"

int commandlist(char *command, FILE *fpin, FILE *fprint)
{ 	
  int nyok=1;
  if      (strcmp(command,"algebra")==0) c_algebra(fpin,fprint); 
  else if (strcmp(command,"alias")==0) c_alias(fpin,fprint); 
  else if (strcmp(command,"areaflowint")==0) c_areaflowint(fpin,fprint); 
  else if (strcmp(command,"arraydelete")==0) c_arraydelete(fpin,fprint); 
  else if (strcmp(command,"arraydump")==0) c_arraydump(fpin,fprint); 
  else if (strcmp(command,"arraydumpmore")==0) c_arraydumpmore(fpin,fprint); 
  else if (strcmp(command,"arrayhowto")==0) c_arrayhowto(fpin,fprint); 
  else if (strcmp(command,"arraylist")==0) c_arraylist(fpin,fprint); 
  else if (strcmp(command,"arrayread")==0) c_arrayread(fpin,fprint); 
  else if (strcmp(command,"aveijk")==0) c_aveijk(fpin,fprint); 
  else if (strcmp(command,"bijrhsmarv")==0) c_bijrhsmarv(fpin,fprint); 
  else if (strcmp(command,"bijrhsmarvs")==0) c_bijrhsmarvs(fpin,fprint); 
  else if (strcmp(command,"bijwallmarv")==0) c_bijwallmarv(fpin,fprint);  
  else if (strcmp(command,"bijwallsplat")==0) c_bijwallsplat(fpin,fprint); 
  else if (strcmp(command,"blocksetabc")==0) c_blocksetabc(fpin,fprint);
  else if (strcmp(command,"coefadd")==0) c_coefadd(fpin,fprint);
  else if (strcmp(command,"coefconv")==0) c_coefconv(fpin,fprint);
  else if (strcmp(command,"coefconvstep")==0) c_coefconvstep(fpin,fprint);
  else if (strcmp(command,"coefcplus")==0) c_coefcplus(fpin,fprint);
  else if (strcmp(command,"coefdt")==0) c_coefdt(fpin,fprint);
  else if (strcmp(command,"coeffix")==0) c_coeffix(fpin,fprint);
  else if (strcmp(command,"coefinit")==0) c_coefinit(fpin,fprint);
  else if (strcmp(command,"coefrhs")==0) c_coefrhs(fpin,fprint);
  else if (strcmp(command,"coefvisc")==0) c_coefvisc(fpin,fprint);
  else if (strcmp(command,"coefviscstep")==0) c_coefviscstep(fpin,fprint);
  else if (strcmp(command,"coefzero")==0) c_coefzero(fpin,fprint);
  else if (strcmp(command,"comment")==0) c_comment(fpin,fprint);
  else if (strcmp(command,"constant")==0) c_constant(fpin,fprint);
  else if (strcmp(command,"contcpcdu")==0) c_contcpcdu(fpin,fprint);
  else if (strcmp(command,"contcpcexit")==0) c_contcpcexit(fpin,fprint);
  else if (strcmp(command,"contcpcfixed")==0) c_contcpcfixed(fpin,fprint);
  else if (strcmp(command,"contcpcinlet")==0) c_contcpcinlet(fpin,fprint);
  else if (strcmp(command,"contdu")==0) c_contdu(fpin,fprint);
  else if (strcmp(command,"contrhsexit")==0) c_contrhsexit(fpin,fprint);
  else if (strcmp(command,"contrhsp")==0) c_contrhsp(fpin,fprint);
  else if (strcmp(command,"contrhsu")==0) c_contrhsu(fpin,fprint);
  else if (strcmp(command,"copy")==0) c_copy(fpin,fprint);
  else if (strcmp(command,"copypart")==0) c_copypart(fpin,fprint);
  else if (strcmp(command,"copy_mtop")==0) c_copy_mtop(fpin,fprint);
  else if (strcmp(command,"copy_ptom")==0) c_copy_ptom(fpin,fprint);
  else if (strcmp(command,"cvdcinit")==0) c_cvdcinit(fpin,fprint);
  else if (strcmp(command,"cvdcparm")==0) c_cvdcparm(fpin,fprint);
  else if (strcmp(command,"cvdcreset")==0) c_cvdcreset(fpin,fprint);
  else if (strcmp(command,"ddtall")==0) c_ddtall(fpin,fprint);
  else if (strcmp(command,"edit")==0) c_edit(fpin,fprint);
  else if (strcmp(command,"editabcd")==0) c_editabcd(fpin,fprint);
  else if (strcmp(command,"eqnsolvebij")==0) c_eqnsolvebij(fpin,fprint);
  else if (strcmp(command,"eqnsolvep")==0) c_eqnsolvep(fpin,fprint);
  else if (strcmp(command,"eqnsolves")==0) c_eqnsolves(fpin,fprint);
  else if (strcmp(command,"eqnupdatem")==0) c_eqnupdatem(fpin,fprint);
  else if (strcmp(command,"eqppts2ppts")==0) c_eqppts2ppts(fpin,fprint);
  else if (strcmp(command,"eqpts2gpts")==0) c_eqpts2gpts(fpin,fprint);
  else if (strcmp(command,"exitinit")==0) c_exitinit(fpin,fprint);
  else if (strcmp(command,"function")==0) c_function(fpin,fprint);
  else if (strcmp(command,"geom8print")==0) c_geom8print(fpin,fprint);
  else if (strcmp(command,"geomcprint")==0) c_geomcprint(fpin,fprint);
  else if (strcmp(command,"gpts2eqpts")==0) c_gpts2eqpts(fpin,fprint);
  else if (strcmp(command,"gradprop")==0) c_gradprop(fpin,fprint);
  else if (strcmp(command,"gridcorner")==0) c_gridcorner(fpin,fprint);
  else if (strcmp(command,"gridfrommefp")==0) c_gridfrommefp(fpin,fprint);
  else if (strcmp(command,"gridmatch")==0) c_gridmatch(fpin,fprint);
  else if (strcmp(command,"gridtomefp")==0) c_gridtomefp(fpin,fprint);
  else if (strcmp(command,"gridvarmod")==0) c_gridvarmod(fpin,fprint);
  else if (strcmp(command,"if")==0) c_if(fpin,fprint);
  else if (strcmp(command,"inletinit")==0) c_inletinit(fpin,fprint);
  else if (strcmp(command,"inletreset")==0) c_inletreset(fpin,fprint);
  else if (strcmp(command,"interp_gtom")==0) c_interp_gtom(fpin,fprint);
  else if (strcmp(command,"interp_gtop")==0) c_interp_gtop(fpin,fprint);
  else if (strcmp(command,"interp_ptod")==0) c_interp_ptod(fpin,fprint);
  else if (strcmp(command,"interp_ptog")==0) c_interp_ptog(fpin,fprint);
  else if (strcmp(command,"lineoutput")==0) c_lineoutput(fpin,fprint);
  else if (strcmp(command,"lineoutputijk")==0) c_lineoutputijk(fpin,fprint);
  else if (strcmp(command,"mirror")==0) c_mirror(fpin,fprint);
  else if (strcmp(command,"momcam")==0) c_momcam(fpin,fprint);
  else if (strcmp(command,"momcamddt")==0) c_momcamddt(fpin,fprint);
  else if (strcmp(command,"momrhsr")==0) c_momrhsr(fpin,fprint);
  else if (strcmp(command,"omrhscoakley")==0) c_omrhscoakley(fpin,fprint);
  else if (strcmp(command,"omrhsmarv")==0) c_omrhsmarv(fpin,fprint);
  else if (strcmp(command,"omwall")==0) c_omwall(fpin,fprint);
  else if (strcmp(command,"omwallcoakley")==0) c_omwallcoakley(fpin,fprint);
  else if (strcmp(command,"omwallmarv")==0) c_omwallmarv(fpin,fprint);
  else if (strcmp(command,"omwallmarvs")==0) c_omwallmarvs(fpin,fprint);
  else if (strcmp(command,"ppreset")==0) c_ppreset(fpin,fprint);
  else if (strcmp(command,"ppts2eqppts")==0) c_ppts2eqppts(fpin,fprint);
  else if (strcmp(command,"prinstress")==0) c_prinstress(fpin,fprint);
  else if (strcmp(command,"print")==0) c_print(fpin,fprint);
  else if (strcmp(command,"prints")==0) c_prints(fpin,fprint);
  else if (strcmp(command,"printscoef")==0) c_printscoef(fpin,fprint);
  else if (strcmp(command,"qturbrhs")==0) c_qturbrhs(fpin,fprint);
  else if (strcmp(command,"rmsminmax")==0) c_rmsminmax(fpin,fprint);
  else if (strcmp(command,"set_contar")==0) c_set_contar(fpin,fprint);
  else if (strcmp(command,"set_cpda")==0) c_set_cpda(fpin,fprint);
  else if (strcmp(command,"set_cpflop")==0) c_set_cpflop(fpin,fprint);
  else if (strcmp(command,"set_cpsleep")==0) c_set_cpsleep(fpin,fprint);
  else if (strcmp(command,"set_gbij")==0) c_set_gbij(fpin,fprint);
  else if (strcmp(command,"set_pkdk")==0) c_set_pkdk(fpin,fprint);
  else if (strcmp(command,"set_srate")==0) c_set_srate(fpin,fprint);
  else if (strcmp(command,"set_volcont")==0) c_set_volcont(fpin,fprint);
  else if (strcmp(command,"set_volmom")==0) c_set_volmom(fpin,fprint);
  else if (strcmp(command,"set_wherefw")==0) c_set_wherefw(fpin,fprint);
  else if (strcmp(command,"set_wherep")==0) c_set_wherep(fpin,fprint);
  else if (strcmp(command,"set_xyzdouble")==0) c_set_xyzdouble(fpin,fprint);
  else if (strcmp(command,"valueat")==0) c_valueat(fpin,fprint);
  else if (strcmp(command,"varinit")==0) c_varinit(fpin,fprint);
  else if (strcmp(command,"varmatch")==0) c_varmatch(fpin,fprint);
  else if (strcmp(command,"varupdate")==0) c_varupdate(fpin,fprint);
  else if (strcmp(command,"vint2eqpts")==0) c_vint2eqpts(fpin,fprint);
  else if (strcmp(command,"visccoakley")==0) c_visccoakley(fpin,fprint);
  else if (strcmp(command,"viscmarv")==0) c_viscmarv(fpin,fprint);
  else if (strcmp(command,"viscmarvheat")==0) c_viscmarvheat(fpin,fprint);
  else if (strcmp(command,"viscqnoise")==0) c_viscqnoise(fpin,fprint);
  else if (strcmp(command,"walldist")==0) c_walldist(fpin,fprint);
  else if (strcmp(command,"wallflux")==0) c_wallflux(fpin,fprint);
  else if (strcmp(command,"wallnorm")==0) c_wallnorm(fpin,fprint);

  /* experimental routines for partial models being tested */
  
  else if (strcmp(command,"bijrhsmarvex")==0) ec_bijrhsmarvex(fpin,fprint); 
  else if (strcmp(command,"qbijles")==0) ec_qbijles(fpin,fprint);
  else if (strcmp(command,"rayleigh")==0) ec_rayleigh(fpin,fprint);
  else if (strcmp(command,"viscles")==0) ec_viscles(fpin,fprint);
  else if (strcmp(command,"w2wlineflat")==0) ec_w2wlineflat(fpin,fprint);
  
/* optional plot commands */  
  
  else if (strcmp(command,"abcmask")==0) pc_abcmask(fpin,fprint);
  else if (strcmp(command,"bar")==0) pc_bar(fpin,fprint);
  else if (strcmp(command,"giflist")==0) pc_giflist(fpin,fprint);
  else if (strcmp(command,"gridinfo")==0) pc_gridinfo(fpin,fprint);
  else if (strcmp(command,"image")==0) pc_image(fpin,fprint);
  else if (strcmp(command,"iplaneint")==0) pc_iplaneint(fpin,fprint);
  else if (strcmp(command,"keyword")==0) pc_keyword(fpin,fprint);
  else if (strcmp(command,"lineplot")==0) pc_lineplot(fpin,fprint);
  else if (strcmp(command,"picture")==0) pc_picture(fpin,fprint);
  else if (strcmp(command,"xytocgrid")==0) pc_xytocgrid(fpin,fprint);
  else if (strcmp(command,"xyzicut")==0) pc_xyzicut(fpin,fprint);
  
  else nyok=0;
  
  return nyok;
}
