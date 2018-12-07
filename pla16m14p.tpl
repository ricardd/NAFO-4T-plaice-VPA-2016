//modified barom2b.tpl in codseal\zap\Mproj\admb
// modified m3use2c - replace initial fixed M with a prior for initial M
// modified codmr1 for plaice
// rvq upper bound on logrvq of 0
// add mobile sentinel, all q's estimated
//2M model with plus group lowered to 16+
DATA_SECTION
  init_int syr;  init_int lyr;  //use actual year 1976 2015
  init_int pyr;
  init_int sage;  init_int lage; // 4 16p
  init_int rvsage; init_int rvlage; // 4 16p
  init_int mssage; init_int mslage; // 4 16p
  init_int rv1_yr1;  init_int rv1_yr2;  // 1976 2002
  init_int rv2_yr1;  init_int rv2_yr2;  // 2004 2015
  init_int ms_yr1;  init_int ms_yr2;  // 2003 2011
  init_number rvmon; 
  init_number msmon; 
  init_matrix rv1dat(rvsage,rvlage,rv1_yr1,rv1_yr2)
  init_matrix rv2dat(rvsage,rvlage,rv2_yr1,rv2_yr2)
  init_matrix msdat(mssage,mslage,ms_yr1,ms_yr2)
  init_matrix C_ay(sage,lage,syr,lyr)
  init_number Minit1;
  init_number Minit2;
  init_matrix wt_a(sage,lage,syr,lyr)
  init_vector wt_term(sage,lage) 
  init_matrix mat_a(sage,lage,syr,lyr)
  init_vector mat_term(sage,lage) 
  init_int Mphz
  init_int phzMdev
  init_number Mdevstd1
  init_number Mdevstd2
  init_number mid_age  //age of split 1
  init_number frdat // F ratio, ages 20p vs age 19, if not estimated
  init_matrix CwtAge(sage,lage,syr,lyr);
  init_vector PRav(sage,lage);
  init_vector Cpro(1,4);
  init_number baranovIter;
  vector age(4,16);
  !! age.fill_seqadd(4,1); 
  init_int icheck
  
  number Cplus2
  
    //projection variables for output
      vector Mpro(sage,lage);    // average of last 5 yr
      vector Spro1(lyr+1,pyr);
      vector Spro2(lyr+1,pyr);
      vector Spro3(lyr+1,pyr);
      vector Spro4(lyr+1,pyr);
      matrix Fpro(1,4,lyr+1,pyr); //  Ftg 
      matrix Npro1(sage,lage,lyr+1,pyr);    
      matrix Npro2(sage,lage,lyr+1,pyr);    
      matrix Npro3(sage,lage,lyr+1,pyr);    
      matrix Npro4(sage,lage,lyr+1,pyr);    

  
 LOC_CALCS
   if(icheck!=12345)
   {
     cout<<" data entry error Mdevstd1 "<<Mdevstd1<<endl;
     cout<<" data entry error Mdevstd2 "<<Mdevstd2<<endl;
     cout<<" data entry error midage "<<mid_age<<endl;
     cout<<" data entry error frdat "<<frdat<<endl;
     cout<<" data entry error CwtAge "<<CwtAge<<endl;
     cout<<" data entry error wt_term "<<wt_term<<endl;
     cout<<" data entry error mat_term "<<mat_term<<endl;
     cout<<" data entry error rvmon"<<rvmon<<endl;
     cout<<" data entry error Cpro "<<Cpro<<endl;
     cout<<" data entry error icheck "<<icheck<<endl;
     exit(1);
   }
 END_CALCS
   


PARAMETER_SECTION
  init_bounded_vector Nterm_y(sage+1,lage,0,15,1)
  init_bounded_vector log_qrv(rvsage,rvlage,-5,0.0,1)
  init_bounded_vector log_qms(mssage,mslage,-5,0.5,1)
  init_bounded_number Mavg1(0.2,1.8,Mphz)
  init_bounded_number Mavg2(0.05,0.8,Mphz)
  init_bounded_vector Mdev1(syr+1,lyr,-5,5,phzMdev) // start in 1977 
  init_bounded_vector Mdev2(syr+1,lyr,-5,5,phzMdev)
  init_bounded_vector log_sig1(rvsage,rvlage,-3.0,1.5,4)
  init_bounded_vector log_sig2(mssage,mslage,-3.0,1.5,4)
//  init_bounded_vector frest(lyr-5,lyr,0.1,2.0,1)
//  init_bounded_vector log_S50(1,2,0.6,2.0,phzS50);
//  init_bounded_vector log_S95step(1,2,0.0,2.0,phzS95);

  
  number Fplus2
//  number qRV;
//  number qMS;
//  vector S50(1,2);
//  vector S95(1,2);
//  matrix selage(1,2,sage,lage);
  
  matrix M_ay(sage,lage,syr,lyr)
  matrix F_ay(sage,lage,syr,lyr)
  matrix N_ay(sage,lage,syr,lyr)
  matrix N_ay_Sept(sage,lage,syr,lyr)
  vector rvq(rvsage,rvlage)  

  matrix yhat1(rvsage,rvlage,rv1_yr1,rv1_yr2)
  matrix RV1resid_ay(rvsage,rvlage,rv1_yr1,rv1_yr2)
  matrix RV1pred_ay(rvsage,rvlage,rv1_yr1,rv1_yr2)

  matrix yhat2(rvsage,rvlage,rv2_yr1,rv2_yr2)
  matrix RV2resid_ay(rvsage,rvlage,rv2_yr1,rv2_yr2)
  matrix RV2pred_ay(rvsage,rvlage,rv2_yr1,rv2_yr2)

  vector msq(mssage,mslage)
  matrix MSyhat(mssage,mslage,ms_yr1,ms_yr2)
  matrix MSresid_ay(mssage,mslage,ms_yr1,ms_yr2)
  matrix MSpred_ay(mssage,mslage,ms_yr1,ms_yr2)

  matrix msdatSept(mssage,mslage,ms_yr1,ms_yr2)

  vector sig1(rvsage,rvlage);
  vector sig2(mssage,mslage);
  vector SSB(syr,lyr+1);
  vector SSBsep(syr,lyr+1);
  matrix SSBatage(sage,lage,syr,lyr+1);
  vector R(syr,lyr);
  vector biom10p(syr,lyr+1);
  vector abun49(syr,lyr);
  vector abun10p(syr,lyr);  
  vector F10p(syr,lyr);
  vector F49(syr,lyr);
  vector Nsep49(rv1_yr1,rv2_yr2);
  vector Nsep10p(rv1_yr1,rv2_yr2);
  number n2term // n2term is average age-4 abundance lyr-2 to lyr 
  sdreport_vector sdbiom10p(syr,lyr+1);
  sdreport_vector sdm1(syr,lyr);
  sdreport_vector sdm2(syr,lyr);
  sdreport_vector sdf1(syr,lyr);
  sdreport_vector sdf2(syr,lyr);
  likeprof_number m176;
  likeprof_number m276;
  vector fcomp(1,20);
  objective_function_value f;


PROCEDURE_SECTION
  //cout<<" start procedure "<<endl;
  fcomp=0.;
  get_m();
  //cout<<" finished get m "<<endl;
  do_dynamics_baranov();
  //cout<<" finished dynamics "<<endl;
  sdbiom10p=biom10p/1000;
  sdm1=M_ay(sage);
  sdm2=M_ay(mid_age);
  sdf1=F49;
  sdf2=F10p;
  m176=Mavg1;
  m276=Mavg2;
  //calc_selage();
  //cout<<" finished sel_age "<<endl;
  get_likelihoods();
  //cout<<" finished likelihoods "<<endl;
  //cout<<" Fmat "<<F_ay<<endl;
   //exit(1);
  f=sum(fcomp);
  if(mceval_phase() )
  {
    count_eval++;
       projection();
       mcout1<<M_ay(sage)<<" "<<M_ay(mid_age)<<endl;
       mcout2<<F49<<" "<<F10p<<endl;
       mcout3<<SSB<<" "<<biom10p<<endl;
       mcout10<<SSBsep<<" "<<SSBsep<<endl;
       mcout8<<SSBatage<<endl;
       mcout9<<N_ay_Sept<<endl;
       mcout6<<R<<endl;
       mcout4<<Nsep49<<endl;
       mcout5<<Nsep10p<<endl;
       mcout7<<rvq<<endl;
       mcoutSp1   << Spro1   << endl;
       mcoutSp2   << Spro2   << endl;
       mcoutSp3   << Spro3   << endl;
       mcoutSp4   << Spro4   << endl;
  }
  compute_N_Sept();

FUNCTION get_m
  int iy,ia;
  M_ay=Mavg1;
  for (ia=mid_age; ia<=lage; ia++) M_ay(ia)=Mavg2;
     
  if(active(Mdev1))
  {
    for (iy=syr+1; iy<=lyr; iy++) 
     {
       for (ia=sage; ia<mid_age; ia++) M_ay(ia,iy)=M_ay(ia,iy-1)*mfexp(Mdev1(iy));
       for (ia=mid_age; ia<=lage; ia++) M_ay(ia,iy)=M_ay(ia,iy-1)*mfexp(Mdev2(iy));
     }
  }


FUNCTION dvariable getF_backward(const dvariable& N,const double& C,const dvariable& M)
  dvariable Ctmp=C;
  dvariable ntmp=N*mfexp(M)+Ctmp*mfexp(M/2.0);
  dvariable fest=log(ntmp/N)-M;
  dvariable ztmp;
  dvariable Pcat;
  dvariable df;
  for (int iter=1; iter<=5; iter++)
  {
    ztmp=fest+M;
    Pcat=fest/ztmp*(mfexp(ztmp)-1)*N;
    df=N/ztmp*(mfexp(ztmp)*(1-fest/ztmp+fest)-1 + fest/ztmp);
    fest+= -1*(Pcat-Ctmp)/df;
  }
  return fest; 


FUNCTION do_dynamics_baranov
  int iy,ia;
  dvariable ztmp;
  SSB=0;
  SSBsep=0;
  SSBatage=0;
  R=0;
  biom10p=0;
  abun49=0;
  abun10p=0;
  F49=0;
  F10p=0;
  n2term=0;
  Nsep49=0.;
  Nsep10p=0.;
  
// step 1: change frdat to frest(lyr) if estimating ratio (3 occurrences)
  for (ia=sage; ia<=lage-2; ia++)
  { 
    F_ay(ia,lyr)=getF_backward(mfexp(Nterm_y(ia+1)),C_ay(ia,lyr),M_ay(ia,lyr));
    N_ay(ia,lyr)=mfexp(Nterm_y(ia+1))*mfexp(M_ay(ia,lyr)+F_ay(ia,lyr));
  }
  Cplus2=C_ay(lage-1,lyr)+C_ay(lage,lyr);
  Fplus2=getF_backward(mfexp(Nterm_y(lage)),Cplus2,M_ay(lage,lyr));
  F_ay(lage-1,lyr)=Fplus2*(C_ay(lage,lyr)+frdat*C_ay(lage-1,lyr))/(frdat*Cplus2);
  ztmp=F_ay(lage-1,lyr)+M_ay(lage-1,lyr);
  N_ay(lage-1,lyr)=(C_ay(lage-1,lyr)*ztmp)/(F_ay(lage-1,lyr)*(1-mfexp(-ztmp)));
  F_ay(lage,lyr)=frdat*F_ay(lage-1,lyr);
  ztmp=F_ay(lage,lyr)+M_ay(lage,lyr);
  N_ay(lage,lyr)=(C_ay(lage,lyr)*ztmp)/(F_ay(lage,lyr)*(1-mfexp(-ztmp)));     
  
// step 2: change frdat to frest(iy) if estimating ratio  (3 occurrences)
  for (iy=lyr-1; iy>=lyr-5; iy--)
  {
    for (ia=sage; ia<=lage-2; ia++) 
    {
      F_ay(ia,iy)=getF_backward(N_ay(ia+1,iy+1),C_ay(ia,iy),M_ay(ia,iy));
      N_ay(ia,iy)= N_ay(ia+1,iy+1)*mfexp(M_ay(ia,iy)+F_ay(ia,iy));
    }
    Cplus2=C_ay(lage-1,iy)+C_ay(lage,iy);
    Fplus2=getF_backward(N_ay(lage,iy+1),Cplus2,M_ay(lage,iy));
    F_ay(lage-1,iy)=Fplus2*(C_ay(lage,iy)+frdat*C_ay(lage-1,iy))/(frdat*Cplus2);
    ztmp=F_ay(lage-1,iy)+M_ay(lage-1,iy);
    N_ay(lage-1,iy)=(C_ay(lage-1,iy)*ztmp)/(F_ay(lage-1,iy)*(1-mfexp(-ztmp)));
    F_ay(lage,iy)=frdat*F_ay(lage-1,iy);
    ztmp=F_ay(lage,iy)+M_ay(lage,iy);
    N_ay(lage,iy)=(C_ay(lage,iy)*ztmp)/(F_ay(lage,iy)*(1-mfexp(-ztmp)));     
  }

// step 3  
  for (iy=lyr-6; iy>=syr; iy--)
  {
    for (ia=sage; ia<=lage-2; ia++) 
    {
      F_ay(ia,iy)=getF_backward(N_ay(ia+1,iy+1),C_ay(ia,iy),M_ay(ia,iy));
      N_ay(ia,iy)= N_ay(ia+1,iy+1)*mfexp(M_ay(ia,iy)+F_ay(ia,iy));
    }
    Cplus2=C_ay(lage-1,iy)+C_ay(lage,iy);
    Fplus2=getF_backward(N_ay(lage,iy+1),Cplus2,M_ay(lage,iy));
    F_ay(lage-1,iy)=Fplus2*(C_ay(lage,iy)+frdat*C_ay(lage-1,iy))/(frdat*Cplus2);
    ztmp=F_ay(lage-1,iy)+M_ay(lage-1,iy);
    N_ay(lage-1,iy)=(C_ay(lage-1,iy)*ztmp)/(F_ay(lage-1,iy)*(1-mfexp(-ztmp)));
    F_ay(lage,iy)=frdat*F_ay(lage-1,iy);
    ztmp=F_ay(lage,iy)+M_ay(lage,iy);
    N_ay(lage,iy)=(C_ay(lage,iy)*ztmp)/(F_ay(lage,iy)*(1-mfexp(-ztmp)));     
  }  
// calculate SSB biom10p F49 F914 F15p abun49 abun914 abun15p vectors  
  for (iy=syr; iy<=lyr; iy++)
  {
    R(iy)=N_ay(sage,iy);
    for (ia=sage; ia<=lage; ia++) SSB(iy)+=N_ay(ia,iy)*wt_a(ia,iy)*mat_a(ia,iy);
    for (ia=sage; ia<=lage; ia++) SSBatage(ia,iy)=N_ay(ia,iy)*wt_a(ia,iy)*mat_a(ia,iy);
    for (ia=10; ia<=lage; ia++) biom10p(iy)+=N_ay(ia,iy)*wt_a(ia,iy);
    for (ia=sage; ia<=9; ia++) abun49(iy)+=N_ay(ia,iy);
    for (ia=10; ia<=lage; ia++) abun10p(iy)+=N_ay(ia,iy);
  }  
  for (iy=syr; iy<=lyr; iy++)
  {
    for (ia=sage; ia<=9; ia++) F49(iy)+=F_ay(ia,iy)*N_ay(ia,iy)/abun49(iy); 
    for (ia=10; ia<=lage; ia++) F10p(iy)+=F_ay(ia,iy)*N_ay(ia,iy)/abun10p(iy); 
  }  
  
  n2term=(N_ay(sage,lyr)+N_ay(sage,lyr-1)+N_ay(sage,lyr-2))/3;
  SSB(lyr+1)=n2term*wt_term(sage)*mat_term(sage);
  for (ia=sage+1; ia<=lage; ia++) SSB(lyr+1)+=mfexp(Nterm_y(ia))*wt_term(ia)*mat_term(ia);
  for (ia=10; ia<=lage; ia++) biom10p(lyr+1)+=mfexp(Nterm_y(ia))*wt_term(ia);

//calculate Nseps and SSBsep
  for (iy=rv1_yr1; iy<=rv2_yr2; iy++)
  {
    for (ia=sage; ia<=lage; ia++) SSBsep(iy)+=N_ay(ia,iy)*wt_a(ia,iy)*mat_a(ia,iy)*mfexp(-.083333*rvmon*(F_ay(ia,iy)+M_ay(ia,iy)));
    for (ia=sage; ia<=9; ia++) Nsep49(iy)+=N_ay(ia,iy)*mfexp(-.083333*rvmon*(F_ay(ia,iy)+M_ay(ia,iy)));
    for (ia=10; ia<=lage; ia++) Nsep10p(iy)+=N_ay(ia,iy)*mfexp(-.083333*rvmon*(F_ay(ia,iy)+M_ay(ia,iy)));
  }


 
FUNCTION get_likelihoods
   int ia,iy;
    yhat1.initialize();
    yhat2.initialize();
    MSyhat.initialize();
    fcomp.initialize();
    dvariable sumSq;
    int nobs;
    

//RV
   for (ia=rvsage; ia<=rvlage; ia++) 
   {
      rvq(ia) = mfexp(log_qrv(ia));
      sig1(ia) = mfexp(log_sig1(ia));
      nobs=0.;
      sumSq=0.;
      for (iy=rv1_yr1; iy<=rv1_yr2; iy++)
      { 
	    yhat1(ia,iy) = N_ay(ia,iy)*mfexp(-.083333*rvmon*(F_ay(ia,iy)+M_ay(ia,iy)));
            RV1pred_ay(ia) = yhat1(ia)*rvq(ia);
	    RV1resid_ay(ia,iy) = log(rv1dat(ia,iy)/RV1pred_ay(ia,iy));
	    sumSq += pow( RV1resid_ay(ia,iy), 2);
	    nobs  += 1;
      }
      fcomp(1) += nobs*log(sig1(ia))+(sumSq/(2.*square(sig1(ia))));
   }
   for (ia=rvsage; ia<=rvlage; ia++) 
   {
      nobs=0.;
      sumSq=0.;
      for (iy=rv2_yr1; iy<=rv2_yr2; iy++)
      { 
	    yhat2(ia,iy) = N_ay(ia,iy)*mfexp(-.083333*rvmon*(F_ay(ia,iy)+M_ay(ia,iy)));
            RV2pred_ay(ia) = yhat2(ia)*rvq(ia);
	    RV2resid_ay(ia,iy) = log(rv2dat(ia,iy)/RV2pred_ay(ia,iy));
	    sumSq += pow( RV2resid_ay(ia,iy), 2);
	    nobs  += 1;
      }
      fcomp(2) += nobs*log(sig1(ia))+(sumSq/(2.*square(sig1(ia))));
   }

//Mobile Sentinel   
   for (ia=mssage; ia<=mslage; ia++)
   {
      nobs=0.;
      sumSq=0.;
      msq(ia) = mfexp(log_qms(ia));
      sig2(ia) = mfexp(log_sig2(ia));
      for (iy=ms_yr1; iy<=ms_yr2; iy++)
      { 
	 MSyhat(ia,iy) = N_ay(ia,iy)*mfexp(-.083333*msmon*(F_ay(ia,iy)+M_ay(ia,iy)));
         MSpred_ay(ia) = MSyhat(ia)*msq(ia);
	 MSresid_ay(ia,iy) = log(msdat(ia,iy)/MSpred_ay(ia,iy));
	 sumSq += pow( MSresid_ay(ia,iy), 2);
	 nobs  += 1;
      }
      fcomp(3) += nobs*log(sig2(ia))+(sumSq/(2.*square(sig2(ia))));
   }

// M random walk terms
   fcomp(4)=0.5*norm2(Mdev1)/square(Mdevstd1);   
   fcomp(5)=0.5*norm2(Mdev2)/square(Mdevstd2); 

// normal prior for Mavg
   fcomp(7)=get_prior_value(1,Mavg1,Minit1,0.05);
   fcomp(8)=get_prior_value(1,Mavg2,Minit2,0.05);

 //  RV_std_resid=sqrt(norm2(RVresid_ay-mean(RVresid_ay))/size_count(RVresid_ay));


FUNCTION dvariable get_prior_value(const int& type,const dvariable& par,const double& mode,const double& std)
// 1 is normal; 2 is lognormal; 3 is robust normal; 4 is robust lognormal; 5 is log-uniform (fixed so not breen)
   dvariable tmp;
   switch (type)
  {
     case 1:  tmp=0.5*square((par-mode)/std)+log(std);  break;
     case 2:  tmp= 0.5*square(log(par/mode)/std + 0.5*std)  + log(par) ; break;
     case 3:  tmp= -log(mfexp(-0.5*square((par-mode)/std))+0.01) + log(std);  break;
     case 4:  tmp= -log(mfexp(-1*square(log(mode/par)+0.5*square(std))
                      /(2*square(std)))+0.01) + log(par) + log(std);       break;
     case 5:  tmp= log(par);  break;  // assumes that log-uniform parameters are parameterized on log scale!!!!!
     default: cout<<"Prior not defined for type "<<type<<endl;      break;
  }
  return tmp;

FUNCTION compute_N_Sept
// loop over year
    for (int y=syr; y<=lyr; y++)
    {
              for (int a=sage; a<=lage; a++)
              {
              N_ay_Sept(a,y) = N_ay(a,y)*mfexp(-((rvmon/12)*F_ay(a,y))-((rvmon/12)*M_ay(a,y)));
              }
    }
// same thing for MS survey
   for (int y=ms_yr1; y<=ms_yr2; y++)
   {
          for(int a=mssage; a<=mslage; a++)
          {
          msdatSept(a,y) = msdat(a,y) * mfexp(-((1/12)*F_ay(a,y))-((1/12)*M_ay(a,y)));
          }
   }




FUNCTION  double getFforward( const int& t, const int& p, const int& nIter, const dvector& Np);
  {
  double fp;
  double Jp;
  double Fest;
  dvector Bp(sage,lage);
  dvector tmp(sage,lage);
  dvector Zap(sage,lage); 
  dvector ZapNew(sage,lage);  
  
  // Initialize Z to current vector of Mta
  ZapNew.initialize(); 
  // Initial approximation of F...
  // Selected biomass
  Bp = elem_prod(elem_prod( Np,column(CwtAge,lyr)), PRav );
  Fest = Cpro(p)/sum( Bp );
  ZapNew = Mpro + PRav*Fest;
  // refine F for fisheries (surveys not critical)
  for( int i=1; i<=nIter; i++ )
  {
    // Total mortality
    Zap=ZapNew; ZapNew=Mpro;
    // Predicted catch given F
    tmp    = elem_div( elem_prod( Bp*Fest,1.-exp(-Zap) ), Zap );
    //cout << "predCatch = " << tmp << endl;
    // Function value: difference of pred - obs catch
    fp =  sum(tmp) - Cpro(p); tmp=0;
    // Jacobian
    dvector tmp1 = elem_div( Bp, Zap );
    dvector tmp2 = elem_prod( PRav, Zap )*Fest;
    dvector tmp3 = 1. - mfexp( -Zap );
    dvector tmp4 = elem_prod( elem_div( PRav*Fest, Zap ), 1.-exp(-Zap)  );
    
    tmp = elem_prod( tmp1, tmp2 + tmp3 - tmp4 );
    Jp = sum(tmp); tmp=0;
    Fest -= fp/Jp;
    ZapNew += PRav*Fest; 
    //cout <<"t = " << t<< " iter = "<< i << " f = "<< f <<" J = "<< J << " f/J = " << f/J << endl;   
    //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, 1) << endl;   
  }
  return Fest;
  }

FUNCTION projection
 {
   int i; int j; int k; int p;
   dmatrix usewate(sage,lage,lyr+1,pyr);
  // randomly pick weight-at-age vectors from last 5 years
   for (j=lyr+1; j<=pyr; j++)
   {
     //k=20.*randu(rng);
     k=5.*randu(rng);
     for (i=sage; i<=lage; i++) usewate(i,j)=wt_a(i,lyr-k);
   }
  
  // assign M to average of last 5 years
  Mpro = 0.0;
  for (i=sage; i<=lage; i++) 
    {
     for (j=lyr-4; j<=lyr; j++) Mpro(i) += value(M_ay(i,j))/5;
    }
    
  // make data variables for required dvariables
  // nT a relic from Sean model = last year
  double rnyr;
  dvector MnT(sage,lage);
  dvector NnT(sage,lage);
  dvector FnT(sage,lage);
  dvector ssb(syr,lyr);
  dvector recruits(syr,lyr);
  dvector recrate(1,10);
  dvector rRateUse(lyr+1,pyr);
  
  MnT = value(column(M_ay,lyr));
  NnT = value(column(N_ay,lyr));
  ssb = value(SSB(syr,lyr));
  FnT = value(column(F_ay,lyr));
  recruits = value(N_ay(sage));
  for (i=1; i<=10; i++) recrate(i) = recruits(lyr-10+i)/ssb(lyr-14+i);
  // get recrates to use
     //rnyr = 9.;
     rnyr = 5.;
     for (j=lyr+1; j<=pyr; j++)
     {
       k = rnyr*randu(rng) + 1;
       rRateUse(j) = recrate(k);
     }
  
//*************************************** Projection 1, Catch = 0 t  ************************************
  // PROJECTION 1, Cpro = 0 t
  Spro1 = 0.0;
  Npro1 = 0.0;
  
  // first projection year
  Npro1(sage,lyr+1) = rRateUse(lyr+1) *ssb(lyr-3);  
  for (i=sage+1; i<=lage-1; i++) Npro1(i,lyr+1) = NnT(i-1)*mfexp(-MnT(i-1)-FnT(i-1));
  Npro1(lage,lyr+1) =  NnT(lage-1)*mfexp(-MnT(lage-1)-FnT(lage-1));
  Npro1(lage,lyr+1) +=  NnT(lage)*mfexp(-MnT(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro1(lyr+1) += Npro1(i,lyr+1)*usewate(i,lyr+1)*mat_a(i,lyr);
  Fpro(1,lyr+1)=getFforward( lyr+1, 1, baranovIter, column(Npro1,lyr+1));
  
  for( j=lyr+2; j<=lyr+4; j++)
  {
    Npro1(sage,j) = rRateUse(j) *ssb(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro1(i,j) = Npro1(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(1,j-1)));
    Npro1(lage,j) =  Npro1(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(1,j-1)));
    Npro1(lage,j) +=  Npro1(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(1,j-1)));
    for (i=sage; i<=lage; i++) Spro1(j) += Npro1(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(1,j)=getFforward( j, 1, baranovIter, column(Npro1,j));
  }
  
  // remaining projection years
  for( j=lyr+5; j<=pyr; j++)
  {
    Npro1(sage,j) = rRateUse(j) *Spro1(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro1(i,j) = Npro1(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(1,j-1)));
    Npro1(lage,j) =  Npro1(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(1,j-1)));
    Npro1(lage,j) +=  Npro1(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(1,j-1)));
    for (i=sage; i<=lage; i++) Spro1(j) += Npro1(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(1,j)=getFforward( j, 1, baranovIter, column(Npro1,j));
  }

 
//*************************************** END PROJECTION1 ************************************

//*************************************** Projection 2, Catch = 100.0 t  ************************************
  // PROJECTION 2, Cpro = 100 t
  Spro2 = 0.0;
  Npro2 = 0.0;
  
  // first projection year
  Npro2(sage,lyr+1) = rRateUse(lyr+1) *ssb(lyr-3);  
  for (i=sage+1; i<=lage-1; i++) Npro2(i,lyr+1) = NnT(i-1)*mfexp(-MnT(i-1)-FnT(i-1));
  Npro2(lage,lyr+1) =  NnT(lage-1)*mfexp(-MnT(lage-1)-FnT(lage-1));
  Npro2(lage,lyr+1) +=  NnT(lage)*mfexp(-MnT(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro2(lyr+1) += Npro2(i,lyr+1)*usewate(i,lyr+1)*mat_a(i,lyr);
  Fpro(2,lyr+1)=getFforward( lyr+1, 2, baranovIter, column(Npro2,lyr+1));
  
  for( j=lyr+2; j<=lyr+4; j++)
  {
    Npro2(sage,j) = rRateUse(j) *ssb(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro2(i,j) = Npro2(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(2,j-1)));
    Npro2(lage,j) =  Npro2(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(2,j-1)));
    Npro2(lage,j) +=  Npro2(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(2,j-1)));
    for (i=sage; i<=lage; i++) Spro2(j) += Npro2(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(2,j)=getFforward( j, 2, baranovIter, column(Npro2,j));
  }
  
  // remaining projection years
  for( j=lyr+5; j<=pyr; j++)
  {
    Npro2(sage,j) = rRateUse(j) *Spro2(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro2(i,j) = Npro2(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(2,j-1)));
    Npro2(lage,j) =  Npro2(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(2,j-1)));
    Npro2(lage,j) +=  Npro2(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(2,j-1)));
    for (i=sage; i<=lage; i++) Spro2(j) += Npro2(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(2,j)=getFforward( j, 2, baranovIter, column(Npro2,j));
  }

 
//*************************************** END PROJECTION2 ************************************

//*************************************** Projection 3, Catch = 250.0 t  ************************************
  // PROJECTION 3, Cpro = 250 t
  Spro3 = 0.0;
  Npro3 = 0.0;
  
  // first projection year
  Npro3(sage,lyr+1) = rRateUse(lyr+1) *ssb(lyr-3);  
  for (i=sage+1; i<=lage-1; i++) Npro3(i,lyr+1) = NnT(i-1)*mfexp(-MnT(i-1)-FnT(i-1));
  Npro3(lage,lyr+1) =  NnT(lage-1)*mfexp(-MnT(lage-1)-FnT(lage-1));
  Npro3(lage,lyr+1) +=  NnT(lage)*mfexp(-MnT(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro3(lyr+1) += Npro3(i,lyr+1)*usewate(i,lyr+1)*mat_a(i,lyr);
  Fpro(3,lyr+1)=getFforward( lyr+1, 3, baranovIter, column(Npro3,lyr+1));
  
  for( j=lyr+2; j<=lyr+4; j++)
  {
    Npro3(sage,j) = rRateUse(j) *ssb(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro3(i,j) = Npro3(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(3,j-1)));
    Npro3(lage,j) =  Npro3(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(3,j-1)));
    Npro3(lage,j) +=  Npro3(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(3,j-1)));
    for (i=sage; i<=lage; i++) Spro3(j) += Npro3(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(3,j)=getFforward( j, 3, baranovIter, column(Npro3,j));
  }
  
  // remaining projection years
  for( j=lyr+5; j<=pyr; j++)
  {
    Npro3(sage,j) = rRateUse(j) *Spro3(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro3(i,j) = Npro3(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(3,j-1)));
    Npro3(lage,j) =  Npro3(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(3,j-1)));
    Npro3(lage,j) +=  Npro3(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(3,j-1)));
    for (i=sage; i<=lage; i++) Spro3(j) += Npro3(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(3,j)=getFforward( j, 3, baranovIter, column(Npro3,j));
  }

 
//*************************************** END PROJECTION3 ************************************

//*************************************** Projection 4, Catch = 0 t, F = half of last 5 years  ************************************
  Mpro = 0.0;
  for (i=sage; i<=lage; i++) 
    {
     for (j=lyr-4; j<=lyr; j++) Mpro(i) += value(M_ay(i,j))/10;
    }
  // PROJECTION 4, Cpro = 250 t
  Spro4 = 0.0;
  Npro4 = 0.0;
  
  // first projection year
  Npro4(sage,lyr+1) = rRateUse(lyr+1) *ssb(lyr-3);  
  for (i=sage+1; i<=lage-1; i++) Npro4(i,lyr+1) = NnT(i-1)*mfexp(-MnT(i-1)-FnT(i-1));
  Npro4(lage,lyr+1) =  NnT(lage-1)*mfexp(-MnT(lage-1)-FnT(lage-1));
  Npro4(lage,lyr+1) +=  NnT(lage)*mfexp(-MnT(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro4(lyr+1) += Npro4(i,lyr+1)*usewate(i,lyr+1)*mat_a(i,lyr);
  Fpro(4,lyr+1)=getFforward( lyr+1, 4, baranovIter, column(Npro4,lyr+1));
  
  for( j=lyr+2; j<=lyr+4; j++)
  {
    Npro4(sage,j) = rRateUse(j) *ssb(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro4(i,j) = Npro4(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(4,j-1)));
    Npro4(lage,j) =  Npro4(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(4,j-1)));
    Npro4(lage,j) +=  Npro4(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(4,j-1)));
    for (i=sage; i<=lage; i++) Spro4(j) += Npro4(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(4,j)=getFforward( j, 4, baranovIter, column(Npro4,j));
  }
  
  // remaining projection years
  for( j=lyr+5; j<=pyr; j++)
  {
    Npro4(sage,j) = rRateUse(j) *Spro4(j-4);  
    for (i=sage+1; i<=lage-1; i++) Npro4(i,j) = Npro4(i-1,j-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(4,j-1)));
    Npro4(lage,j) =  Npro4(lage-1,j-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(4,j-1)));
    Npro4(lage,j) +=  Npro4(lage,j-1)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(4,j-1)));
    for (i=sage; i<=lage; i++) Spro4(j) += Npro4(i,j)*usewate(i,j)*mat_a(i,lyr);
    Fpro(4,j)=getFforward( j, 4, baranovIter, column(Npro4,j));
  }

 
//*************************************** END PROJECTION4 ************************************

 }





GLOBALS_SECTION
  #include <admodel.h>
  ofstream mcout1("mcoutsM.tmp");     //M4-9 M10+
  ofstream mcout2("mcoutsF.tmp");     //F4-9 1014 15+ 
  ofstream mcout3("mcoutsSSB.tmp");   
  ofstream mcout8("mcoutsSSBatage.tmp");   
  ofstream mcout10("mcoutsSSBsep.tmp");   
  ofstream mcout4("mcoutNs1.tmp");
  ofstream mcout5("mcoutNs2.tmp");
  ofstream mcout6("mcoutR.tmp");
  ofstream mcout7("mcoutRVq.tmp");
  ofstream mcoutSp1("mcoutSp1.dat");
  ofstream mcoutSp2("mcoutSp2.dat");
  ofstream mcoutSp3("mcoutSp3.dat");
  ofstream mcoutSp4("mcoutSp4.dat");
  ofstream mcout9("mcoutsNatageSept.tmp");   

  int count_eval;
  int iseed=1337;
  random_number_generator rng(iseed);

TOP_OF_MAIN_SECTION
  arrmblsize=1200000;
  int ndepvar=400;
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(ndepvar);
  

REPORT_SECTION
 int ia,iy;
   
  report<<"#popnJan"<<endl<<N_ay<<endl;
  report<<"#popnSept"<<endl<<N_ay_Sept<<endl;
  report<<"#F"<<endl<<F_ay<<endl;
  report<<"#gulfrvdat1"<<endl<<rv1dat<<endl;
  report<<"#gulfrvdat2"<<endl<<rv2dat<<endl; 
  report<<"#RV1res"<<endl<<RV1resid_ay<<endl;
  report<<"#RV2res"<<endl<<RV2resid_ay<<endl;
  report<<"#MSres"<<endl<<MSresid_ay<<endl;
  report<<"#vPar2010"<<endl<<Nterm_y<<endl;
  report<<"#vSSB"<<endl<<SSB<<endl;
  report<<"#SSBatage"<<endl<<SSBatage<<endl;
  report<<"#biom10p"<<endl<<biom10p<<endl;
  report<<"#abun49"<<endl<<abun49<<endl;
  report<<"#abun10p"<<endl<<abun10p<<endl;
  report<<"#F49"<<endl<<F49<<endl;
  report<<"#F10p"<<endl<<F10p<<endl;
  report<<"#nsep1"<<endl<<Nsep49<<endl;
  report<<"#nsep2"<<endl<<Nsep10p<<endl;
  report<<"#m1"<<endl<<M_ay(sage)<<endl;
  report<<"#m2"<<endl<<M_ay(mid_age)<<endl;
  report<<"#rvq"<<endl<<rvq<<endl; 
  report<<"#msq"<<endl<<msq<<endl;
  report<<"#msdatSept"<<endl<<msdatSept<<endl;
  report<<"#RVsig"<<endl<<sig1<<endl; 
  report<<"#MSsig"<<endl<<sig2<<endl; 
  report<<"#fcomp"<<endl<<fcomp<<endl;
  report<<"#vSSBsep"<<endl<<SSBsep<<endl;
  report<<"#maxGrad"<<endl<<objective_function_value::gmax<<endl; 
//  report<<"#Qrv"<<endl<<mfexp(log_qrv)<<endl; 
//  report<<"#Qms"<<endl<<mfexp(log_qms)<<endl; 
 
  
// ofstream out("runt",ios::app);
// if(last_phase())  out<<"M1term "<<M1term<<" M2term "<<M2term<<" "<<Bfinal<<" "<<Ffinal59<<" fcomp "<<fcomp<<" Tot "<<f<<endl;

