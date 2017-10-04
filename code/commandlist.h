int commandlist(char *command, FILE *fpin, FILE *fprint);

/* define all commands as void */
/* commands for experimental test models listed below */
/* followed by plot package commands */

void c_algebra(FILE *fpin, FILE *fprint);
void c_alias(FILE *fpin, FILE *fprint); /* in arrays.c */
void c_areaflowint(FILE *fpin, FILE *fprint);
void c_arraydelete(FILE *fpin, FILE *fprint);
void c_arraydump(FILE *fpin, FILE *fprint);
void c_arraydumpmore(FILE *fpin, FILE *fprint); /* in c_arraydump */ 
void c_arrayhowto(FILE *fpin, FILE *fprint);
void c_arraylist(FILE *fpin, FILE *fprint); /* in arrays.c */
void c_arrayread(FILE *fpin, FILE *fprint);
void c_aveijk(FILE *fpin, FILE *fprint);
void c_bijrhsmarv(FILE *fpin, FILE *fprint);
void c_bijrhsmarvex(FILE *fpin, FILE *fprint);
void c_bijrhsmarvs(FILE *fpin, FILE *fprint); /* in c_bijrhsmarv */
void c_bijwallmarv(FILE *fpin, FILE *fprint);
void c_bijwallsplat(FILE *fpin, FILE *fprint);
void c_blocksetabc(FILE *fpin, FILE *fprint);  /* in blocksubs */
void c_coefadd(FILE *fpin, FILE *fprint);
void c_coefconv(FILE *fpin, FILE *fprint);
void c_coefconvstep(FILE *fpin, FILE *fprint);  /* in c_coefconv.c */
void c_coefcplus(FILE *fpin, FILE *fprint);
void c_coefdt(FILE *fpin, FILE *fprint);
void c_coeffix(FILE *fpin, FILE *fprint);
void c_coefinit(FILE *fpin, FILE *fprint);
void c_coefrhs(FILE *fpin, FILE *fprint);
void c_coefvisc(FILE *fpin, FILE *fprint);
void c_coefviscstep(FILE *fpin, FILE *fprint); /* in c_coefvisc.c */
void c_coefzero(FILE *fpin, FILE *fprint);
void c_comment(FILE *fpin, FILE *fprint);
void c_constant(FILE *fpin, FILE *fprint);
void c_contcpcdu(FILE *fpin, FILE *fprint);
void c_contcpcexit(FILE *fpin, FILE *fprint); /* in exitsubs.c */
void c_contcpcfixed(FILE *fpin, FILE *fprint);
void c_contcpcinlet(FILE *fpin, FILE *fprint); /* in inletsubs.c */
void c_contdu(FILE *fpin, FILE *fprint);
void c_contrhsexit(FILE *fpin, FILE *fprint); /* in exitsubs.c */
void c_contrhsp(FILE *fpin, FILE *fprint);
void c_contrhsu(FILE *fpin, FILE *fprint);
void c_copy(FILE *fpin, FILE *fprint);
void c_copypart(FILE *fpin, FILE *fprint);
void c_copy_mtop(FILE *fpin, FILE *fprint);
void c_copy_ptom(FILE *fpin, FILE *fprint);
void c_cvdcinit(FILE *fpin, FILE *fprint); /* in cvdcsubs.c */
void c_cvdcparm(FILE *fpin, FILE *fprint); /* in cvdcsubs.c */
void c_cvdcreset(FILE *fpin, FILE *fprint); /* in cvdcsubs.c */
void c_ddtall(FILE *fpin, FILE *fprint);
void c_edit(FILE *fpin, FILE *fprint);
void c_editabcd(FILE *fpin, FILE *fprint);
void c_eqnsolvebij(FILE *fpin, FILE *fprint);
void c_eqnsolvep(FILE *fpin, FILE *fprint);
void c_eqnsolves(FILE *fpin, FILE *fprint);
void c_eqnupdatem(FILE *fpin, FILE *fprint);
void c_eqppts2ppts(FILE *fpin, FILE *fprint);
void c_eqpts2gpts(FILE *fpin, FILE *fprint);
void c_exitinit(FILE *fpin, FILE *fprint); /* in exitsubs.c */
void c_function(FILE *fpin, FILE *fprint);
void c_geom8print(FILE *fpin, FILE *fprint); /* in geom8subs.c */
void c_geomcprint(FILE *fpin, FILE *fprint); /* in geomcsubs.c */
void c_gpts2eqpts(FILE *fpin, FILE *fprint);
void c_gradprop(FILE *fpin, FILE *fprint);
void c_gridcorner(FILE *fpin, FILE *fprint);
void c_gridfrommefp(FILE *fpin, FILE *fprint);
void c_gridmatch(FILE *fpin, FILE *fprint);
void c_gridtomefp(FILE *fpin, FILE *fprint);
void c_gridvarmod(FILE *fpin, FILE *fprint);
void c_if(FILE *fpin, FILE *fprint);
void c_inletinit(FILE *fpin, FILE *fprint); /* in inletsubs.c */
void c_inletreset(FILE *fpin, FILE *fprint); /* in inletsubs.c */
void c_interp_gtom(FILE *fpin, FILE *fprint);
void c_interp_gtop(FILE *fpin, FILE *fprint);
void c_interp_ptod(FILE *fpin, FILE *fprint);
void c_interp_ptog(FILE *fpin, FILE *fprint);
void c_lineoutput(FILE *fpin, FILE *fprint);
void c_lineoutputijk(FILE *fpin, FILE *fprint);
void c_mirror(FILE *fpin, FILE *fprint);
void c_momcam(FILE *fpin, FILE *fprint);
void c_momcamddt(FILE *fpin, FILE *fprint);
void c_momrhsr(FILE *fpin, FILE *fprint);
void c_omrhscoakley(FILE *fpin, FILE *fprint);
void c_omrhsmarv(FILE *fpin, FILE *fprint);
void c_omwall(FILE *fpin, FILE *fprint); /* in omwallsubs */
void c_omwallcoakley(FILE *fpin, FILE *fprint); /* in omwallsubs */
void c_omwallmarv(FILE *fpin, FILE *fprint);  /* in omwallsubs */
void c_omwallmarvs(FILE *fpin, FILE *fprint);  /* in omwallsubs */
void c_ppreset(FILE *fpin, FILE *fprint);
void c_ppts2eqppts(FILE *fpin, FILE *fprint);
void c_prinstress(FILE *fpin, FILE *fprint);
void c_print(FILE *fpin, FILE *fprint);
void c_prints(FILE *fpin, FILE *fprint);
void c_printscoef(FILE *fpin, FILE *fprint);
void c_qturbrhs(FILE *fpin, FILE *fprint);
void c_rmsminmax(FILE *fpin, FILE *fprint);
void c_set_contar(FILE *fpin, FILE *fprint); /* in c_set_cpda */
void c_set_cpda(FILE *fpin, FILE *fprint);
void c_set_cpflop(FILE *fpin, FILE *fprint);
void c_set_cpsleep(FILE *fpin, FILE *fprint);
void c_set_gbij(FILE *fpin, FILE *fprint);
void c_set_pkdk(FILE *fpin, FILE *fprint);
void c_set_srate(FILE *fpin, FILE *fprint);
void c_set_volcont(FILE *fpin, FILE *fprint);
void c_set_volmom(FILE *fpin, FILE *fprint);
void c_set_wherefw(FILE *fpin, FILE *fprint);
void c_set_wherep(FILE *fpin, FILE *fprint);
void c_set_xyzdouble(FILE *fpin, FILE *fprint);
void c_valueat(FILE *fpin, FILE *fprint);
void c_varinit(FILE *fpin, FILE *fprint);
void c_varmatch(FILE *fpin, FILE *fprint);
void c_varupdate(FILE *fpin, FILE *fprint);
void c_vint2eqpts(FILE *fpin, FILE *fprint);
void c_visccoakley(FILE *fpin, FILE *fprint);
void c_viscmarv(FILE *fpin, FILE *fprint);
void c_viscmarvheat(FILE *fpin, FILE *fprint);
void c_viscqnoise(FILE *fpin, FILE *fprint);
void c_walldist(FILE *fpin, FILE *fprint);
void c_wallflux(FILE *fpin, FILE *fprint);
void c_wallnorm(FILE *fpin, FILE *fprint);

/* experimental routines for partial models being tested */

void ec_bijrhsmarvex(FILE *fpin, FILE *fprint);
void ec_qbijles(FILE *fpin, FILE *fprint);
void ec_rayleigh(FILE *fpin, FILE *fprint);
void ec_viscles(FILE *fpin, FILE *fprint);
void ec_w2wlineflat(FILE *fpin, FILE *fprint);

/* optional plot commands */  

void pc_abcmask(FILE *fpin, FILE *fprint);
void pc_bar(FILE *fpin, FILE *fprint);
void pc_giflist(FILE *fpin, FILE *fprint);
void pc_gridinfo(FILE *fpin, FILE *fprint);
void pc_image(FILE *fpin, FILE *fprint);
void pc_keyword(FILE *fpin, FILE *fprint);
void pc_lineplot(FILE *fpin, FILE *fprint);
void pc_picture(FILE *fpin, FILE *fprint);
void pc_xyzicut(FILE *fpin, FILE *fprint);
void pc_iplaneint(FILE *fpin, FILE *fprint);
void pc_xytocgrid(FILE *fpin, FILE *fprint);
