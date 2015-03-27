DATA_SECTION

	init_int seed;
	int syr;
	int eyr;
	int nage;
	int nlen;
	int lstp;
	
	LOC_CALCS
		syr = 1;
		eyr = 38;
		nage = 12;
		nlen = 77;
		lstp = 1;
	END_CALCS
	
	init_vector ft(syr,eyr);
	number sigR;
	number tau;
	number m;
	number ahat;
	number ghat;
	number Linf;
	number k;
	number to;
	number cvl;
	number alw;
	number blw;
	number reca;
	number recb;
	number Eo;
	number Ro;
	number reck;
	number Am1;
	number phie;
	number ssb0;
	number test;
	number q;

	vector wt(syr,eyr);
	vector eps(syr,eyr);
// 	vector ft(syr,eyr); 
	vector age(1,nage);
	vector vbt(syr,eyr);
// 	vector sbt(syr,eyr
	vector ct(syr,eyr);
	vector bt(syr,eyr);
	
	vector len(1,nlen);
	vector va(1,nage);
	vector vul(1,nage);
	vector lxo(1,nage);
	vector wa(1,nage);
	vector fec(1,nage);
	vector la(1,nage);
	vector Sa(1,nage);
	vector iwt(syr,eyr-1);
	vector iyr(1,eyr)
	vector z1(1,nlen);
	vector z2(1,nlen);
	vector std(1,nage);
	

	matrix P_la(1,nage,1,nlen);
	matrix P_al(1,nlen,1,nage);
	matrix Nat(syr,eyr+1,1,nage);
	matrix zt(syr,eyr,1,nage);
	matrix cat(syr,eyr,1,nage);
	matrix pcat(syr,eyr,1,nage);
	matrix cal(syr,eyr,1,nlen);
	
	LOC_CALCS
		sigR = 0.2;
		tau = 0.1;
		ahat = 3.643583;
		ghat = 1.143614	;
		m = 0.23;
		to = 0;
		cvl = 0.09;
		alw = 0.0000238;
		blw = 2.7671;
		Ro = 2000;
		Linf = 76.464;
		k = 0.09;
		reck = 7.111111;
		random_number_generator rng(seed);
		wt.fill_randn(rng);
		wt*=sigR;
		eps.fill_randn(rng);
		eps*=tau;
		age.fill_seqadd(1,1);
		len.fill_seqadd(lstp,lstp);
		Am1=nage-1;
		q=1;
		iwt =1;
		iyr.fill_seqadd(1,1);
		
	END_CALCS
	
PARAMETER_SECTION
	objective_function_value no_f; 
	  
 	LOC_CALCS
		incidence_functions();
		output_data();
		output_true();
		exit(1);
	END_CALCS

PROCEDURE_SECTION

FUNCTION incidence_functions
	
	//  biological processes that will not be directly estimated in the assessment
	z1(1,nlen); 
	z2(1,nlen);
	std(1,nage);

	Sa = exp(-m);

// 	la.initialize();
	la = Linf*(1.-exp(-k*(age-to))); //average length at age

	std = la * cvl; 

	lxo(1) = 1.; //first age
	
	for(int a=2;a<=nage;a++)
	{
		lxo(a) = lxo(1)*pow(Sa(a-1),age(a-1)); // proportion of individuals at age surviving M only
	}
	
	lxo(nage) /= 1.-Sa(Am1); // age plus group
	
	wa = alw * pow(la,blw); //weight at age
	fec = elem_prod(wa,plogis(age,2.448848,0.568183)); // wight at age-weight at 50% maturity
 	va = plogis(age,ahat,ghat);

	phie = sum(elem_prod(lxo,fec)); 
	
	reca = reck/phie; 
	recb = (reck - 1.)/(Ro*phie); 
	ssb0 = Ro*phie;
	test = (reca*ssb0)/(1.+recb*ssb0);
	
// 	cout << test << endl;
// 	exit(1);
	
	
 	for( int a = 1; a <= nage; a++ )
	{
		z1 = (( len - lstp * 0.5 )-la( a ))/std( a );//calculate for each length bin the probs tof being of specific ages
		z2 = (( len + lstp * 0.5 )-la( a ))/std( a );
		for( int b=1; b<= nlen; b++ )
		{
		P_la( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b )); // calculates the proportion of being of a given age given your length
		}
	}
	
	P_al = trans(P_la); //transpose matrix to length by age

	
// 	  for( int i=1;i<=nage; i++ ) 
// 	    Nat(syr,i) = Ro*exp(-(i-1)*m) * exp(wt(i));
// nor now	
	    Nat(syr) = Ro*lxo;
	  
//  	  ft = 0.0001;
	  for(int i=syr;i<=eyr;i++)
	  {
	    
	    zt(i) = m+ft(i)*va;
	    double sbt = fec * Nat(i);

	    Nat(i+1,1) = (reca*sbt)/(1.+recb*sbt)*exp(wt(i));
	    Nat(i+1)(2,nage) = ++elem_prod(Nat(i)(1,nage-1),exp(-zt(i)(1,nage-1)));
	    Nat(i+1,nage) += Nat(i,nage)*exp(-zt(i,nage));

	    dvector t1 = elem_div(ft(i)*va,zt(i));
	    dvector t2 = elem_prod( 1.-exp(-zt(i)), Nat(i) ) ;
	    cat(i) = elem_prod(t1,t2); // check this
	    pcat(i) = cat(i)/sum(cat(i));
	    // true catch at length
// 	    cal(i) = P_al*pcat(i); //   P_al(77,12) X pcat(38*12)
	    cal(i) = P_al*cat(i); //   P_al(77,12) X pcat(38*12)
	  }
	  
	    for(int i=syr;i<=eyr;i++)
	  {
	  vbt = q * Nat.sub(syr,eyr)*elem_prod(wa,va) * exp(eps(i)); // cpue
	  bt = Nat.sub(syr,eyr)* wa * exp(eps(i)); // survey
	   }

	  ct = cat*wa;
// 	 cout << bt<< endl;
// 	 exit(1);

	  
FUNCTION output_data
	ofstream ofs("perujmsra.dat");
	ofs<<"# syr " << endl << syr <<endl;
	ofs<<"# eyr " << endl << eyr <<endl;
	ofs<<"# nage "<< endl << nage <<endl;
	ofs<<"# nlen "<< endl << nlen <<endl;
	ofs<<"# lstp "<< endl << lstp <<endl;
	ofs<<"# SR function " << endl << 1 <<endl;
	ofs<<"# m " << endl << m <<endl;
	ofs<<"# alw " << endl << alw <<endl;
	ofs<<"# blw "<< endl << blw <<endl;
	ofs<<"# wmat "<< endl << 1.240343  <<endl;
	ofs<<"# mat50  "<< endl << 2.448848 <<endl;
	ofs<<"# #matsd " << endl << 0.568183 <<endl;
	ofs<<"# ahat " << endl << ahat <<endl;
	ofs<<"# ghat "<< endl << ghat <<endl;
	ofs<<"# vul "<< endl << va <<endl;
	ofs<<"# nyt "<< endl << eyr <<endl;
	ofs<<"# iyr " << endl << iyr <<endl;
	ofs<<"# yt " << endl << bt <<endl;
	ofs<<"# cal "<< endl << cal <<endl;
	ofs<<"# ilinf "<< endl << Linf <<endl;
	ofs<<"# ik "<< endl << k <<endl;
	ofs<<"# it0 " << endl << to <<endl;
	ofs<<"# icvl " << endl << cvl <<endl;
	ofs<<"# ireck "<< endl << reck <<endl;
	ofs<<"# iRo "<< endl << Ro <<endl;
	ofs<<"# iwt "<< endl << iwt <<endl;
	ofs<<"# cv_it " << endl << tau <<endl;
	ofs<<"# sigR " << endl << sigR <<endl;
	ofs<<"# sigVul " << endl << 0.2 <<endl;
	ofs<<"# phz_reck "<< endl << -2 <<endl;
	ofs<<"# phz_growth  "<< endl << -4  <<endl;
	ofs<<"# use_prior  "<< endl << 1 <<endl;
	ofs<<"# dend " << endl << 999 <<endl;
	
	
FUNCTION output_true
	  
	ofstream ofs("true_data_lsra.rep");
	ofs<<"true_ct" << endl << ct <<endl;
	ofs<<"true_ut" << endl << ft <<endl;
	ofs<<"true_Nat" << endl << Nat.sub(syr,eyr) <<endl;
	ofs<<"true_cal" << endl << cal.sub(syr,eyr) <<endl;