	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	//#include<MyLikelihoods.cpp>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <perujmsra.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  syr.allocate("syr");
  eyr.allocate("eyr");
  nage.allocate("nage");
  nlen.allocate("nlen");
  lstp.allocate("lstp");
  SR.allocate("SR");
  m.allocate("m");
  alw.allocate("alw");
  blw.allocate("blw");
  wmat.allocate("wmat");
  mat50.allocate("mat50");
  matsd.allocate("matsd");
  ahat.allocate("ahat");
  ghat.allocate("ghat");
  va.allocate(1,nage,"va");
  nyt.allocate("nyt");
  iyr.allocate(1,nyt,"iyr");
  survB.allocate(1,nyt,"survB");
  Clt.allocate(syr,eyr,1,nlen,"Clt");
  ilinf.allocate("ilinf");
  ik.allocate("ik");
  to.allocate("to");
  icvl.allocate("icvl");
  ireck.allocate("ireck");
  iRo.allocate("iRo");
  iwt.allocate(syr,eyr-1,"iwt");
  cv_it.allocate("cv_it");
  sigR.allocate("sigR");
  sigVul.allocate("sigVul");
  phz_reck.allocate("phz_reck");
  phz_growth.allocate("phz_growth");
  use_prior.allocate("use_prior");
  dend.allocate("dend");
		
		if( dend != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
  age.allocate(1,nage);
  len.allocate(1,nlen);
		age.fill_seqadd( 1, 1 );
		len.fill_seqadd( lstp, lstp );
		Am1=nage-1;
		tiny=1.e-20;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_Linf.allocate(phz_growth,"log_Linf");
  log_k.allocate(phz_growth,"log_k");
  log_cvl.allocate(phz_growth,"log_cvl");
  log_reck.allocate(phz_reck,"log_reck");
  log_Ro.allocate(1,"log_Ro");
 log_Linf=log(ilinf);
 log_k=log(ik);
 log_cvl=log(icvl);
 log_reck=log(ireck);
 log_Ro=log(iRo);
  log_wt.allocate(syr,eyr-1,-10.,10.,2,"log_wt");
 log_wt = log(iwt);
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  cvl.allocate("cvl");
  #ifndef NO_AD_INITIALIZE
  cvl.initialize();
  #endif
  reck.allocate("reck");
  #ifndef NO_AD_INITIALIZE
  reck.initialize();
  #endif
  Ro.allocate("Ro");
  #ifndef NO_AD_INITIALIZE
  Ro.initialize();
  #endif
  ssvul.allocate("ssvul");
  #ifndef NO_AD_INITIALIZE
  ssvul.initialize();
  #endif
  zstat.allocate(1,nyt,"zstat");
  #ifndef NO_AD_INITIALIZE
    zstat.initialize();
  #endif
  wt.allocate(syr,eyr-1,"wt");
  #ifndef NO_AD_INITIALIZE
    wt.initialize();
  #endif
  Sa.allocate(1,nage,"Sa");
  #ifndef NO_AD_INITIALIZE
    Sa.initialize();
  #endif
  vul.allocate(1,nage,"vul");
  #ifndef NO_AD_INITIALIZE
    vul.initialize();
  #endif
  Eo.allocate("Eo");
  #ifndef NO_AD_INITIALIZE
  Eo.initialize();
  #endif
  reca.allocate("reca");
  #ifndef NO_AD_INITIALIZE
  reca.initialize();
  #endif
  recb.allocate("recb");
  #ifndef NO_AD_INITIALIZE
  recb.initialize();
  #endif
  la.allocate(1,nage,"la");
  #ifndef NO_AD_INITIALIZE
    la.initialize();
  #endif
  wa.allocate(1,nage,"wa");
  #ifndef NO_AD_INITIALIZE
    wa.initialize();
  #endif
  lxo.allocate(1,nage,"lxo");
  #ifndef NO_AD_INITIALIZE
    lxo.initialize();
  #endif
  fec.allocate(1,nage,"fec");
  #ifndef NO_AD_INITIALIZE
    fec.initialize();
  #endif
  P_la.allocate(1,nage,1,nlen,"P_la");
  #ifndef NO_AD_INITIALIZE
    P_la.initialize();
  #endif
  P_al.allocate(1,nlen,1,nage,"P_al");
  #ifndef NO_AD_INITIALIZE
    P_al.initialize();
  #endif
  Nat.allocate(syr,eyr,1,nage,"Nat");
  #ifndef NO_AD_INITIALIZE
    Nat.initialize();
  #endif
  Ulength.allocate(syr,eyr,1,nlen,"Ulength");
  #ifndef NO_AD_INITIALIZE
    Ulength.initialize();
  #endif
  Uage.allocate(syr,eyr,1,nage,"Uage");
  #ifndef NO_AD_INITIALIZE
    Uage.initialize();
  #endif
  Nlt.allocate(syr,eyr,1,nlen,"Nlt");
  #ifndef NO_AD_INITIALIZE
    Nlt.initialize();
  #endif
  maxUy.allocate(syr,eyr,"maxUy");
  #ifndef NO_AD_INITIALIZE
    maxUy.initialize();
  #endif
  muUl.allocate(1,nlen,"muUl");
  #ifndef NO_AD_INITIALIZE
    muUl.initialize();
  #endif
  psurvB.allocate(syr,eyr,"psurvB");
  #ifndef NO_AD_INITIALIZE
    psurvB.initialize();
  #endif
  q.allocate("q");
  #ifndef NO_AD_INITIALIZE
  q.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
}

void model_parameters::userfunction(void)
{
  nll =0.0;
        trans_parms();
	incidence_functions();
        initialization();
	SRA();
	observation_model();
	objective_function();
}

void model_parameters::trans_parms(void)
{
	Linf = mfexp( log_Linf );
	k = mfexp( log_k );
	cvl = mfexp( log_cvl );
	reck = mfexp( log_reck ) + 1.001;
	Ro = mfexp( log_Ro );
	wt = mfexp( log_wt );
}

void model_parameters::incidence_functions(void)
{
	dvector z1( 1, nlen );
	dvector z2( 1, nlen );
	double zn;
	dvar_vector std( 1, nage );
	la = Linf * ( 1. - mfexp( -k * ( age - to )));
	std = la * cvl;
	Sa = exp(-m);
	lxo( 1 ) = 1.;
	for( int a = 2; a <= nage; a++ )
		lxo( a ) = lxo( 1 ) * pow( Sa( a - 1 ), age( a - 1 ));
		lxo( nage ) /= 1. - Sa( Am1 );
	wa = alw * pow( la, blw );
	fec=elem_prod(wa,plogis(age,mat50,matsd));
	vul = plogis(age,ahat,ghat);
 	for( int a = 1; a <= nage; a++ )
	{
		if( fec( a ) < 0. )	fec( a ) = 0.;
		z1 = (( len - lstp * 0.5 )-value( la( a )))/value( std( a ));
		z2 = (( len + lstp * 0.5 )-value( la( a )))/value( std( a ));
		for( int b=1; b<= nlen; b++ )
		P_la( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b ));
	}
	P_al = trans( P_la );
}

void model_parameters::initialization(void)
{
	Eo = Ro * sum( elem_prod( lxo, fec ));
	reca = reck * Ro / Eo;
	if( SR==1 )		recb =( reck - 1. ) / Eo;
	else			recb = log( reck ) / Eo;
}

void model_parameters::SRA(void)
{
	dvar_matrix Upen( syr, eyr, 1, nlen );
	// INITIAL YEAR
	for( int a = 1; a <= nage; a++ )
		Nat( syr, a ) = Ro * pow( Sa( a ), age( a ) - 1. );	// initial age-structure
		Nat( syr, nage ) /= 1. - Sa( nage );
	for( int b = 1; b <= nlen; b++ )
	{
		Nlt( syr, b ) = Nat( syr ) * P_al( b );			// length-structure in year-1
		Ulength( syr, b ) = Clt( syr, b ) / posfun( Nlt( syr, b ), Clt( syr, b ), fpen);	// exploitation by length
	}
	maxUy( syr ) = max( Ulength( syr ));
	for( int a = 1; a <= nage; a++ )
		Uage( syr, a ) = Ulength( syr ) * P_la( a );		// exploitation by age
	// SUBSEQUENT YEARS
	for( int y = syr + 1; y <= eyr; y++ )
	{
	dvariable sbt = fec * Nat(y - 1);				// eggs in year-y
	if( SR ==1 )
	Nat( y, 1 ) = reca * sbt / ( 1. + recb * sbt) * wt( y - 1 );	// recruitment
	else
	Nat( y, 1 ) = reca * sbt * mfexp( - recb * sbt) * wt( y - 1 );
	// age-distribution post-recruitment
	Nat( y )( 2, nage ) =++ posfun( elem_prod( elem_prod( Nat( y - 1 )( 1, Am1 ), Sa( 1, Am1 )), 1. - Uage( y - 1 )( 1, Am1 )), tiny, fpen );
	Nat( y, nage ) = posfun( Nat( y - 1, nage -1 ) * Sa( nage ) * ( 1. - Uage( y - 1, nage -1)), tiny, fpen); //  += CHECK THE PLUS GROUP
	Nat( y, nage ) += posfun( Nat( y - 1, nage ) * Sa( nage ) * ( 1. - Uage( y - 1, nage )), tiny, fpen); //  += CHECK THE PLUS GROUP
	for( int b = 1; b <= nlen; b++ )
	{
	Nlt( y, b ) = Nat( y ) * P_al( b );						// length-distribution
	Ulength( y, b ) = Clt( y, b ) / posfun( Nlt( y, b ), Clt( y, b ), fpen);	// exploitation by length
	}
	for( int a = 1; a <= nage; a++ )
	Uage( y, a ) = Ulength( y ) * P_la( a );			// exploitation by age
	maxUy( y ) = max( Ulength( y ));				// max exploitation across lengths
	}
	for( int b = 1; b <= nlen; b++ )
	{
	muUl(b) = sum(elem_div( column(Ulength,b),maxUy ) ) / size_count(maxUy);		// mean exploitation rate across years
	}
	for( int y = syr; y <= eyr; y++ )
	Upen( y ) = pow( Ulength( y ) / posfun( maxUy( y ), tiny, fpen ) - muUl, 2. );		// penalty against dramatic changes in vulnerability
 	ssvul = sum( Upen );									// vulnerability penalty
	psurvB = Nat * elem_prod(wa,va);
	dvar_matrix temp_vlt =  elem_div( Clt,Nlt );
}

void model_parameters::observation_model(void)
{
	zstat=log(survB)-log(psurvB(iyr));		// posterior probability distribution of q
	q=mfexp(mean(zstat));
 	zstat -= mean(zstat);				// z-statistic used for calculating MLE of q
}

void model_parameters::objective_function(void)
{
	dvar_vector lvec(1,2);
	lvec.initialize();
	lvec(1)=dnorm(zstat,cv_it);
	lvec(2)=dnorm(wt,sigR);
	dvar_vector pvec(1,4);
	pvec.initialize();
	if(active(log_reck))
	{  
 	dvariable h=reck/(4.+reck);	
	pvec(1)=dbeta((h-0.2)/0.8,1.,1.);
   	}
	if(last_phase())
	{
	pvec(2)=dnorm(log_wt,2.); 	// estimate recruitment deviations with dnorm function
	}
	else
	{
	pvec(3)=100.*norm2(log_wt); 	// assume deviations are at their normal mle analytical estimator  (obs-pred)^2
	}
	pvec(4) = ssvul/sigVul;
	nll = sum(lvec) + sum(pvec)*use_prior;
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	report<<"Ro\n"<<Ro<<endl;
	report<<"bo\n"<<Eo<<endl;
	report<<"kappa\n"<<reck<<endl;
	dvariable h = value(reck/(4.+reck));
	REPORT(h);
	REPORT(Linf);
	REPORT(k);
	REPORT(cvl);
 	REPORT(wt);
	REPORT(reca);
	REPORT(recb);
	REPORT(psurvB);
	REPORT(survB);
	report<<"yt\n"<<survB<<endl;
	report<<"bt\n"<<Nat*wa<<endl;
	dvar_vector wl = alw * pow( len, blw );
	REPORT(wl);
	REPORT(muUl);
	REPORT(maxUy);
	report<<"N\n"<<Nat<<endl;
	REPORT(Ulength);
	REPORT(Uage);
	REPORT(Clt);
	REPORT(Nlt);
 	dvector sbt=value(Nat.sub(syr,eyr)*fec);
 	REPORT(sbt);
	double depletion = sbt(eyr)/sbt(syr);
	REPORT(depletion);
 	ivector yr(syr,eyr);
	yr.fill_seqadd(syr,1); 
	report<<"yr\n"<<yr<<endl;
	report<<"iyr\n"<<iyr<<endl;
	REPORT(age);
	REPORT(len);
	REPORT(q);
	REPORT(Sa);
	REPORT(vul);
	REPORT(la);
	REPORT(va);
	REPORT(wa);
}

void model_parameters::final_calcs()
{
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
		cout<<"*******************************************"<<endl;
		cout<<"--Start time: "<<ctime(&start)<<endl;
		cout<<"--Finish time: "<<ctime(&finish)<<endl;
		cout<<"--Runtime: ";
		cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
		cout<<"*******************************************"<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
