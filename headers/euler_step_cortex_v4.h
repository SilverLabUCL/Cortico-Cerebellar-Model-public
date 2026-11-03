void euler_step_cortex(double t, double hcx[], netpara &para) {

	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);

	// RNN firing rates
	for(int i= 0; i < para.Ncx; i++) {
		para.nuCx[i]=(*para.Ptr_ctxIO)(hcx[i],1.0);
	}

	//	RNN feedback input
	if(para.closedloop==0){
		// open loop: drive network with target signal
		feedback=func_target(t, para);
	}
	else {
		// closed loop: drive network with output signal
		feedback=0;
		for(int i= 0; i < para.Ncx; i++) {
			feedback+=para.nuCx[i]*para.Wout[i];
		}
	}

	// add RNN interference signal
	if(para.intrfflag==1){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
		for(int i= 0; i < para.Ncx; i++) {
			for (int j=0; j<para.Ncomp_intf; j++){
				interference[i]+=input_intf[j]*para.Wfb_interf.at(i,j);
			}
		}
	}

	// recurrent inputs
	for(int i= 0; i < para.Ncx; i++) {
		recurrent[i]=0;
		for (int j=0; j<para.Wcx[i].Nsyn; j++) recurrent[i]+=para.nuCx[para.Wcx[i].idx[j]]*para.Wcx[i].weight[j];
	}

	for(int i= 0; i < para.Ncx; i++) {
		// hcx[i]= (1-dt/tau)*hcx[i]+dt/tau*(recurrent+feedback+interference)+ para.sqDT*sigmaN*xi;
		hcx[i]= para.IIdt_distb[i]*hcx[i]+para.dttau_distrb[i]*(recurrent[i]+feedback*para.Wfb[i]+interference[i]);
	}

	// add noise only during test phase
	if(para.testphaseflag==1){
		for(int i= 0; i < para.Ncx; i++) {
			if(para.noiseflag==0){
				// constant noise on inputs
				hcx[i]= hcx[i]+para.sqDT[i] *gsl_ran_gaussian_ziggurat(para.r,para.sigmaN);
			}
			else if(para.noiseflag==1){
				// noise scaled with sqrt(input)
				hcx[i]= hcx[i]+para.sqDT[i] *gsl_ran_gaussian_ziggurat(para.r,para.sigmaN*sqrt(std::fabs(hcx[i])));
			}
		}
	}

	return;
}
