//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brett van Poorten
//Date:	June 21, 2013; Update: 8 july 2014 
//Purpose: length-based SRA / VPA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	init_int syr;				// start year of survey data
	init_int eyr;				// end year of survey data
	init_int nage;				// number of age-classes
	init_int nlen;				// number of length-bins (assume start at 0)
	init_int lstp;				// width of length-bins
	init_int SR;				// stock-recruit relationship: 1==Beverton-Holt; 2==Ricker
	init_number m;			// natural mortality
	init_number alw;			// multiplier for length-weight relationship
	init_number blw;		// multiplier for length-weight relationship
	init_number wmat;			// initial weight at maturity
	init_number mat50;		// maturity parameter
	init_number matsd;		// maturity parameter
	init_number ahat;		// vulnerability parameter
	init_number ghat;		// vulnerability parameter
	init_vector va(1,nage);			// survey vulnerability
	init_int nyt
	init_ivector iyr(1,nyt); // iyr(syr,nyr)
	init_vector survB(1,nyt);  // yt(syr,nyr)
// 	init_vector survB(syr,eyr);		// survey biomass
	init_matrix Clt(syr,eyr,1,nlen);	// catch at length and year
// !! 	cout<<iyr<<"\n"<<survB<<endl;
// !! 	exit(1);
	

	// what are these for?
	init_number ilinf;
	init_number ik
	init_number to
	init_number icvl;
	init_number ireck
	init_number iRo
 	init_vector iwt(syr,eyr-1);

	init_number cv_it;		// CV for cpue
	init_number sigR;		// sigma R
	init_number sigVul;
	init_int phz_reck;		// phaze for reck
	init_int phz_growth;		// phaze for growth parameters
	init_int use_prior;		// add priors to the likehoods ? (1 or 0)

	init_int dend;
	
	
	LOCAL_CALCS
		
		if( dend != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}

	END_CALCS
	
	vector age(1,nage);			// ages
	vector len(1,nlen);			// middle length of each length bin
	int Am1;				// maximum age minus 1
	number tiny;
	
	LOC_CALCS
		age.fill_seqadd( 1, 1 );
		len.fill_seqadd( lstp, lstp );
		Am1=nage-1;
		tiny=1.e-20;
	END_CALCS

//	!! cout<<"ilinf\t"<<ilinf<<"\n ik\t"<<ik<<"\n icvl\t"<<icvl<<"\n iRo\t"<<iRo<<"\n ireck\t"<<ireck<<"\n lstp\t"<<lstp<<"\n Sa\t"<<Sa<<endl;
//	!! exit(1);

PARAMETER_SECTION
	init_number log_Linf(phz_growth);
	init_number log_k(phz_growth);
	init_number log_cvl(phz_growth);	
//	init_bounded_number log_Linf(1.09,2.48,-1)
//	init_bounded_number log_k(-4.6,-0.69,-1)
// 	init_bounded_number log_cvl(-4.60,-0.91,-3);
	init_number log_reck(phz_reck);
	init_number log_Ro(1);
	
	!! log_Linf=log(ilinf);
	!! log_k=log(ik);
	!! log_cvl=log(icvl);
	!! log_reck=log(ireck);
	!! log_Ro=log(iRo);

	init_bounded_dev_vector log_wt(syr,eyr-1,-10.,10.,2);

 	!! log_wt = log(iwt);


	objective_function_value nll;
	
	number fpen;
	number Linf;					// von Bertalanffy asymptotic length
	number k;					// von Bertalanffy metabolic parameter
	number cvl;					// coefficient of variation in length at age
	number reck;					// Goodyear recruitment compensation parameter
	number Ro;					// unfished recruitment
	number ssvul;
// 	vector zstat(syr,eyr);				// posterior probability distribution of q
	vector zstat(1,nyt);	
	vector wt(syr,eyr-1);				// recruitment anomolies
	vector Sa(1,nage);		// survival-at-age (assume constant across ages)
 	vector vul(1,nage);		// age-specific vulnerabilities
	number Eo;					// unfished egg deposition
	number reca;					// alpha of stock-recruit relationship
	number recb;					// beta of stock-recruit relationship
	vector la(1,nage);				// length-at-age
	vector wa(1,nage);				// weight-at-age
	vector lxo(1,nage);				// unfished survivorship at age
	vector fec(1,nage);				// age-specific fecundity - used for stock-recruit function
	matrix P_la(1,nage,1,nlen);			// probability of being length at age
	matrix P_al(1,nlen,1,nage);			// transpose of above
	matrix Nat(syr,eyr,1,nage);
	matrix Ulength(syr,eyr,1,nlen);
	matrix Uage(syr,eyr,1,nage);
	matrix Nlt(syr,eyr,1,nlen);
	vector maxUy(syr,eyr);
	vector muUl(1,nlen);
 	vector psurvB(syr,eyr);
	number q;

PRELIMINARY_CALCS_SECTION


PROCEDURE_SECTION
        trans_parms();
	incidence_functions();
        initialization();
	SRA();
	observation_model();
	objective_function();

FUNCTION trans_parms
	Linf = mfexp( log_Linf );
	k = mfexp( log_k );
	cvl = mfexp( log_cvl );
	reck = mfexp( log_reck ) + 1.001;
	Ro = mfexp( log_Ro );
	wt = mfexp( log_wt );
	
// 	cout<<"Linf\t"<<Linf<<"\n k\t"<<k<<"\n cvl\t"<<cvl<<"\n reck\t"<<reck<<"\n Ro\t"<<Ro<<"\n wt\n"<<wt<< endl;
// 	exit(1);

FUNCTION incidence_functions
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
	
//	wa = alw * pow( la, 3. );
//	fec = wa - wmat;

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
	
//  	cout<<"la\t"<<la<<"\n wa\t"<<wa<<"\n lxo\t"<<lxo<<"\n fec\t"<<fec<<"\n len\t"<<len<<"\n P_la\n"<<P_la<<"\n P_al\n"<<P_al<<endl;
// 	exit(1);
	
FUNCTION initialization
	Eo = Ro * sum( elem_prod( lxo, fec ));
	reca = reck * Ro / Eo;
	if( SR==1 )		recb =( reck - 1. ) / Eo;
	else			recb = log( reck ) / Eo;

FUNCTION SRA
	dvar_matrix Upen( syr, eyr, 1, nlen );
// 	dvar_matrix Nlt( syr, eyr, 1, nlen );
	
	// INITIAL YEAR
	for( int a = 1; a <= nage; a++ )
		Nat( syr, a ) = Ro * pow( Sa( a ), age( a ) - 1. );	// initial age-structure
		Nat( syr, nage ) /= 1. - Sa( nage );
	
	
	for( int b = 1; b <= nlen; b++ )
	{
		Nlt( syr, b ) = Nat( syr ) * P_al( b );			// length-structure in year-1
		Ulength( syr, b ) = Clt( syr, b ) / posfun( Nlt( syr, b ), Clt( syr, b ), fpen);	// exploitation by length
// 		cout<<"by year\t"<<Clt(syr,b)<<"\t"<<Nlt(syr,b)<<"\t"<<Ulength(syr,b)<<endl;
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
// 	maxUy( y ) = mean( Ulength( y ));	// max exploitation across lengths (mean also works)
	
	}

	for( int b = 1; b <= nlen; b++ )
	{
	muUl(b) = sum(elem_div( column(Ulength,b),maxUy ) ) / size_count(maxUy);		// mean exploitation rate across years
	}
	
	
	for( int y = syr; y <= eyr; y++ )
	Upen( y ) = pow( Ulength( y ) / posfun( maxUy( y ), tiny, fpen ) - muUl, 2. );		// penalty against dramatic changes in vulnerability
 	ssvul = sum( Upen );									// vulnerability penalty
	
	psurvB = Nat * elem_prod(wa,va);
	
// 	Report stuff
	dvar_matrix temp_vlt =  elem_div( Clt,Nlt );


// 	cout<<"Nat\n"<<Nat<<"\n Nlt\n"<<Nlt<<"\n Ulength\n"<<Ulength<<"\n Uage\n"<<Uage<<"\n Upen\n"<<Upen<<endl;
// 	cout<<"muUl\n"<<muUl<<"\n maxUy\n"<<maxUy<<"\n ssvul\n"<<ssvul<<"\n P_al\n"<<P_al<<endl;
// 	exit(1);

FUNCTION observation_model
// 	cout<<psurvB(iyr)<<"\n"<<psurvB<<"\n"<<survB <<endl;
// 	exit(1);
	zstat=log(survB)-log(psurvB(iyr));		// posterior probability distribution of q
//	zstat = log( elem_div( survB, psurvB ));
	q=mfexp(mean(zstat));
 	zstat -= mean(zstat);				// z-statistic used for calculating MLE of q
 
	
FUNCTION objective_function 

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


REPORT_SECTION
	
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
// 	REPORT(Nat);
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

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
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

FINAL_SECTION
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
