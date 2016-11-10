
double beta_from_pz(double zlens, string pzfile)
{
   // shamelessly copied from nicaea's nofz.c
  
   FILE *F;
   char c, cpre;
   int i, n, nread;
   double *par_nz, dummy;
   char str[1024], dummy2[128];

   /* Number of lines */
   F = fopen(pzfile.c_str(), "r");
   if(!F) {
    cerr << "Could not open redshift info file " << pzfile << endl;
    return -1;
   }
   n = 0;
   do {
      c = fgetc(F);
      if (c=='#' && cpre=='\n') n--;
      if (c=='\n') n++;
      cpre = c;
   } while (c!=EOF);
   n--;   /* Last line contains z_n 0 */

   rewind(F);
   
   /* Header line */
   nread = fscanf(F, "%s %s\n", dummy2, str);
   if(strcmp(dummy2, "#")!=0) {
     cerr << "Wrong format of redshift info file " << pzfile << ": First line has to be '# <snofz_t>', e.g. '# hist'" << endl;
     return -1.;  
   }
   
   if(strcmp(str,"single")==0) {
     double zsource;
     
     if(fscanf(F, "%lf\n", &zsource)!=1) {
       cerr << "Error while reading n(z) file: two identical doubles expected (first line)" << endl;
       return -1.;
     }
     cerr << "# single source redshift " << zsource << endl;
     assert(zsource>zlens);
     return angularDiameterDistance(zlens,zsource,1000)/angularDiameterDistance(0,zsource,1000);
   }
   
   else if(strcmp(str,"hist")==0) {
     double z0,z1;  
     if(fscanf(F, "%lf %lf\n", &z0, &z1)!=2) {
       cerr << "Error while reading n(z) file: two doubles expected (first line)" << endl;
       return -1.; 
     }
     cout << "### " << z0 << " " << z1 << endl;
     double zprev=0.;
     double pprev=0.;
     
     double sp=0.;
     double sbeta=0.;
     
     for (i=1; i<n; i++) {
        double znow,pnow;
	
        if(fscanf(F, "%lf %lf\n", &znow, &pnow)!=2) {  /* z_i  n_i */
	  cerr << "Error while reading n(z) file: two doubles expected" << endl;
	  return -1;
	}
	cout << znow << " " << pnow << endl;
	assert(znow>=zprev);
	
	if(zprev+znow>2.*zlens) {
	  sbeta += pprev*angularDiameterDistance(zlens,(zprev+znow)/2.,1000)/angularDiameterDistance(0,(zprev+znow)/2.,1000);
	}
	sp += pprev;
	
	zprev=znow; pprev=pnow;
     }
     if(zprev+z1>2.*zlens) {
	  sbeta += pprev*angularDiameterDistance(zlens,(zprev+z1)/2.,1000)/angularDiameterDistance(0,(zprev+z1)/2.,1000);
     }
     sp += pprev;
     cerr << "mean beta=" << sbeta/sp << endl;
     return sbeta/sp;
   }
   
  return 1.; 
}
