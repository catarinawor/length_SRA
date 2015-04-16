//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Roberto Licandeo and Catraina Wor
//Date:	June 21, 2013; Update: 8 july 2014 
//Purpose: length-based SRA / VPA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	
	// Read new seed in
	int seed;

	LOC_CALCS
		ifstream ifs( "seed.txt" ); // if this file is available
		ifs>>seed; //read in the seed
		seed += 10; // add 10 to the seed?
		ofstream ofs( "seed.txt" ); //put out to seed.txt
		ofs<<seed<<endl; //the new value of the seed
	END_CALCS

	//model dimensions
	init_int syr;
	init_int eyr;
	init_int nage;
	init_int nlen;
	init_int lstp;
	init_int niyr;

	init_int proc_err;
	init_int obs_err;
	
	//true parameter values
	init_number sigR; 	    		// standard deviation for recruitment deviations
	init_number tau;				// standard deviation for observation error
	init_number m;					// natural mortality
	init_number ahat;				// selectivity inflexion point
	init_number ghat;				// selectivity sd
	init_number feca;				// maturity curve inflexion point
	init_number fecg;				// maturity sd
	init_number Linf;				// vb Linfinity
	init_number k;					// vb k
	init_number to;					// vb to
	init_number cvl;				// coefficient of variantion for age at length curve
	init_number alw;				// parameter for length weight relationship
	init_number blw;				// parameter for length weight relationship
	init_number reck;				// recruitment compensarion ratio
	init_number Ro;					// Average unfished recruitment
	init_number q;					// catchability coefficient




	init_vector ft(syr,eyr);		// Fishing mortality pattern
	init_vector iyr(1,niyr);     	// years for which a survey is available

	init_int dend;					// end of file number

		LOCAL_CALCS
		if( dend != 999 )
		{
			cout<<"Error reading true values.\n Fix it."<<endl;
			cout<< ft <<endl;
			cout<< dend <<endl;
			ad_exit(1);
		}
	END_CALCS

	// derived quantities
	number reca;					// BH parameter a
	number recb;					// BH parameter b
	number Eo;						// Average unfished egg production	
	number Am1;						// penultimate age
	number phie;					// Equilibrium fished fecundity
	number sbo;						// Unfished stock spawning biomass
	//number test;					// test to see if S-R function is calculated correctly, should return Ro 
	
	vector wt(syr,eyr);				// Recruitment deviations
	vector eps(syr,eyr); 			// Observation errors deviations
	vector age(1,nage); 			// age vector
	vector vbt(syr,eyr); 			// vulnerable biomass (assuming a single gear)
	vector ct(syr,eyr);  			// catches
	vector bt(syr,eyr);	  			// total biomass
	vector sbt(syr,eyr); 			// spawning biomass
	vector depl(syr,eyr);			// depletion based on spawning biomass
	
	vector len(1,nlen); 			// length classes
	vector va(1,nage); 				// fisheries vulnerability at age
	vector lxo(1,nage); 			// Equilibrium unfished numbers ata age
	vector wa(1,nage); 				// Weight at age
	vector fec(1,nage); 			// fecundity at age
	vector la(1,nage); 				// length at age
	vector Sa(1,nage); 				// Survival at age	
	

	matrix P_la(1,nage,1,nlen); 	// matrix of proportion of age at length
	matrix P_al(1,nlen,1,nage); 	// matrix of proportion of age at length (transpose)
	matrix Nat(syr,eyr+1,1,nage); 	// Numbers at age
	matrix zt(syr,eyr,1,nage); 		// Total mortality
	matrix cat(syr,eyr,1,nage); 	// Catch at age (in numbers)
	matrix pcat(syr,eyr,1,nage); 	// Proportions of catch at age
	matrix cal(syr,eyr,1,nlen); 	// Catch at Length
	
	LOC_CALCS
		random_number_generator rng(seed);
		wt.fill_randn(rng);
		wt*=sigR;
		eps.fill_randn(rng);
		eps*=tau;
		age.fill_seqadd(1,1);
		len.fill_seqadd(lstp,lstp);
		Am1=nage-1;	
	END_CALCS
	
PARAMETER_SECTION
	objective_function_value no_f; 

PRELIMINARY_CALCS_SECTION	  
 	
	incidence_functions();
	output_data();
	output_true();
	exit(1);
	
PROCEDURE_SECTION

FUNCTION incidence_functions
	
	//  biological processes that will not be directly estimated in the assessment
	dvector z1(1,nlen); 				// intermediate steps for calculating proportion of age at length
	dvector z2(1,nlen); 				// intermediate steps for calculating proportion of age at length
	dvector std(1,nage); 					// std for length at age curve

	Sa = exp(-m);

 	la.initialize();
	la = Linf*(1.-exp(-k*(age-to)));  //average length at age
	std = la*cvl;  		  //std for length at age

	lxo(1) = 1.; //first age	
	for(int a = 2 ; a <= nage ; a++)
	{
		lxo(a) = lxo(a-1)*Sa(a-1); // proportion of individuals at age surviving M only
	}
	lxo(nage) /= 1.-Sa(nage); // age plus group
	
	wa = alw * pow(la,blw); //weight at age
	fec = elem_prod(wa,plogis(age,feca,fecg)); // wight at age-weight at 50% maturity
 	va = plogis(age,ahat,ghat);

	phie = lxo*fec; 
	
	reca = reck/phie; 
	recb = (reck - 1.)/(Ro*phie); 
	sbo  = Ro*phie;
	//test = (reca*sbo)/(1.+recb*sbo);

	// 	cout << "Ro is:" << Ro <<endl;	
	// 	cout << "Ro_test is:"<< test <<endl;
	// 	exit(1);
	
	// Calculate proportion of legth at age class

 	for( int a = 1; a <= nage; a++ )
	{
		z1 = (( len - lstp * 0.5 )-la( a ))/std( a );
		z2 = (( len + lstp * 0.5 )-la( a ))/std( a );
		for( int b=1; b<= nlen; b++ )
		{
		P_la( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b )); // calculates the proportion of a given age given your length
		}
	}
	
	P_al = trans(P_la); //transpose matrix to length by age

	Nat(syr) = Ro*lxo;
	  
	for(int i=syr;i<=eyr;i++)
	{	    
	    zt(i) = m+ft(i)*va;

	    sbt(i) = fec * Nat(i);

	    Nat(i+1,1) = (reca*sbt(i))/(1.+recb*sbt(i))*exp(wt(i)*obs_err);
	    Nat(i+1)(2,nage) = ++elem_prod(Nat(i)(1,nage-1),exp(-zt(i)(1,nage-1)));
	    //Nat(i+1,nage) += Nat(i,nage)*exp(-zt(i,nage));
		Nat(i+1,nage) /= 1.- exp(-zt(i,nage));


	    dvector t1 = elem_div(ft(i)*va,zt(i));
	    dvector t2 = elem_prod( 1.-exp(-zt(i)), Nat(i) ) ;
	    cat(i) = elem_prod(t1,t2); 
	    pcat(i) = cat(i)/sum(cat(i));

	    // true catch at length
	    cal(i) = P_al*cat(i); //   P_al(77,12) X pcat(38*12)
	
		vbt(i) = q * Nat(i)*elem_prod(wa,va) * exp(eps(i)*proc_err); // cpue
		bt(i) = Nat(i)* wa * exp(eps(i)*proc_err); 				     // survey
		depl(i) = sbt(i)/sbt(1);
	}

	  ct = cat*wa;

	  
FUNCTION output_data
	dvector iwt(syr,eyr-1);				// Recruitment deviations
	iwt.initialize();

	ofstream ofs("jmsra.dat");
	ofs<<"# syr " << endl << syr <<endl;
	ofs<<"# eyr " << endl << eyr <<endl;
	ofs<<"# nage "<< endl << nage <<endl;
	ofs<<"# nlen "<< endl << nlen <<endl;
	ofs<<"# lstp "<< endl << lstp <<endl;
	ofs<<"# SR function " << endl << 1 <<endl;
	ofs<<"# m " << endl << m <<endl;
	ofs<<"# alw " << endl << alw <<endl;
	ofs<<"# blw "<< endl << blw <<endl;
	ofs<<"# mat50  "<< endl << feca <<endl;
	ofs<<"# matsd " << endl << fecg <<endl;
	ofs<<"# ahat " << endl << ahat <<endl;
	ofs<<"# ghat "<< endl << ghat <<endl;
	ofs<<"# vul "<< endl << va <<endl;
	ofs<<"# nyt "<< endl << eyr <<endl;
	ofs<<"# iyr " << endl << iyr <<endl;
	ofs<<"# yt " << endl << vbt <<endl;
	ofs<<"# cal "<< endl << cal <<endl;
	ofs<<"# ilinf "<< endl << Linf <<endl;
	ofs<<"# ik "<< endl << k <<endl;
	ofs<<"# it0 " << endl << to <<endl;
	ofs<<"# icvl " << endl << cvl <<endl;
	ofs<<"# ireck "<< endl << reck <<endl;
	ofs<<"# iRo "<< endl << Ro <<endl;
	ofs<<"# wt "<< endl << exp(iwt) <<endl;
	ofs<<"# cv_it " << endl << tau <<endl;
	ofs<<"# sigR " << endl << sigR <<endl;
	ofs<<"# sigVul " << endl << 0.4 <<endl;
	ofs<<"# phz_reck "<< endl << 2 <<endl;
	ofs<<"# phz_growth  "<< endl << -4  <<endl;
	ofs<<"# use_prior  "<< endl << 0 <<endl;
	ofs<<"# dend " << endl << 999 <<endl;
	
	
FUNCTION output_true
	  
	ofstream ofs("true_data_lsra.rep");
	ofs<<"true_ct" << endl << ct <<endl;
	ofs<<"true_ut" << endl << ft <<endl;
	ofs<<"true_Nat" << endl << Nat.sub(syr,eyr) <<endl;
	ofs<<"true_cal" << endl << cal.sub(syr,eyr) <<endl;
	ofs<<"true_Ro" << endl << Ro <<endl;
	ofs<<"true_reck" << endl << reck <<endl;
	ofs<<"true_sbt" << endl << sbt <<endl;
	ofs<<"true_depl" << endl << depl <<endl;
	ofs<<"true_q" << endl << q <<endl;	



	





