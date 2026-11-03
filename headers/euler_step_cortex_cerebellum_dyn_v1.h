void euler_step_cortex_cerebellum_dyn(double t, double X[], netpara &para) {

	double GoC;
	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec GCinputs(para.Ngc);
	arma::vec MF(para.Nmf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);

	// RNN firing rates
	for(int i= 0; i < para.Ncx; i++) {
		para.nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// add GC interference signal
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
	}

	// GC rates computed from X
	for (int i=0; i<para.Ngc; i++) {
		para.nuGC[i]=para.gain_vector[i]*(*para.Ptr_gcIO) (X[i+para.Ncx]-para.threshold_vector[i],para.sp_alpha);

		// add GC interference signal
		if(para.intrfflag==2){
			for (int j=0; j<para.Ncomp_intf; j++)	para.nuGC[i]=para.nuGC[i]+input_intf[j]*para.Intf_weights.at(i,j);
		}
	}

	// GC filtering
	if(para.filterflag>0)	para.nuGC=filter_GCs(para.nuGC, para);

	//	RNN feedback input
	if(para.closedloop==0){
		// open loop: drive network with target signal
		feedback=func_target(t, para);
	}
	else {
		// closed loop: drive network with output signal
		feedback=0;
		for(int i= 0; i < para.Ngc; i++) {
			feedback+=para.nuGC[i]*para.Wout_PC[i];
		}
	}

	// MFs = sample from cortical activity
	for (int j=0; j<para.Nmf; j++) 	MF[j]=para.nuCx[para.Ncx_to_MF_idx[j]];

	// Golgi cell
	if(para.gocmodeflag==1 || para.gocmodeflag==2){
		//~ GoC=0;
		//~ for (int i=0; i<para.Ngc; i++) {
			//~ GoC+=para.nuGC[i];
		//~ }
		//~ GoC=GoC/double(para.Ngc);
		GoC=arma::mean(para.nuGC);
	}
	else if(para.gocmodeflag==3){
		GoC=arma::stddev(MF);
	}

	// update GC inputs
	for (int i=0; i<para.Ngc; i++) {
		//	pre-synaptic inputs of GC i
		GCinputs[i]=0;
		for (int j=0; j<para.Nsyn; j++) {
			GCinputs[i]+=MF[para.Wmf[i].idx[j]];
			//~ GCinputs[i]+=Wmf[i].weight[j]*MF[Wmf[i].idx[j]];
		}
		GCinputs[i]=GCinputs[i]/para.sqNsyn;
		if(para.gocmodeflag>0) GCinputs[i]+=-para.Wgoc_distrb[i]*GoC;
	}
	for(int i= 0; i < para.Ngc; i++) {
		X[i+para.Ncx]= para.IIdtGC_distb[i]*X[i+para.Ncx]+para.dttauGC_distb[i]*GCinputs[i];
	}

	// add noise
	if(para.noiseflag==0 && para.testphaseflag==1 && para.sigmaNgc>0){
		// constant noise
		for (int i=0; i<para.Ngc; i++) {
			X[i+para.Ncx]=X[i+para.Ncx] + para.sqDTGC[i]*gsl_ran_gaussian_ziggurat(para.r,para.sigmaNgc);
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

	// update recurrent inputs
	for(int i= 0; i < para.Ncx; i++) {
		recurrent[i]=0;
		for (int j=0; j<para.Wcx[i].Nsyn; j++) recurrent[i]+=para.nuCx[para.Wcx[i].idx[j]]*para.Wcx[i].weight[j];
	}

	for(int i= 0; i < para.Ncx; i++) {
		// X[i]= (1-dt/tau)*X[i]+dt/tau*(...);
		X[i]= para.IIdt_distb[i]*X[i]+para.dttau_distrb[i]*(recurrent[i]+feedback*para.Wfb[i]+interference[i]);
	}

	// add noise only during test phase
	if(para.testphaseflag==1){
		for(int i= 0; i < para.Ncx; i++) {
			if(para.noiseflag==0){
				// constant noise on inputs
				X[i]= X[i]+para.sqDT[i] *gsl_ran_gaussian_ziggurat(para.r,para.sigmaN);
			}
			else if(para.noiseflag==1){
				// noise scaled with sqrt(input)
				X[i]= X[i]+para.sqDT[i] *gsl_ran_gaussian_ziggurat(para.r,para.sigmaN*sqrt(std::fabs(X[i])));
			}
		}
	}


	return;
}
