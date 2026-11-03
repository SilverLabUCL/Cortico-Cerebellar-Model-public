int RK_cortex_cereb_rhs(double t, const double hcx[], double dhcx_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	double GoC;
	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference= arma::zeros<arma::vec>(para.Ncx);
	arma::vec MF(para.Nmf);
	arma::vec nuCx(para.Ncx);
	arma::vec nuGC(para.Ngc);
	arma::vec hGC(para.Ngc);

	// RNN firing rates
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(hcx[i],1.0);
	}

	// MFs = sample from cortical activity
	for (int j=0; j<para.Nmf; j++) 	MF[j]=nuCx[para.Ncx_to_MF_idx[j]];

	// Golgi cell
	if(para.gocmodeflag==1 || para.gocmodeflag==2){
		GoC=arma::mean(nuGC);
	}
	else if(para.gocmodeflag==3){
		GoC=arma::stddev(MF);
	}

	// GC interference signal
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
	}

	// Instantaneous GC inputs
	for (int i=0; i<para.Ngc; i++) {
		//	pre-synaptic inputs of GC i
		hGC[i]=0;
		for (int j=0; j<para.Nsyn; j++) {
			hGC[i]+=MF[para.Wmf[i].idx[j]];
			//~ hGC[i]+=Wmf[i].weight[j]*MF[Wmf[i].idx[j]];
		}
		hGC[i]=hGC[i]/para.sqNsyn;
		if(para.gocmodeflag>0) hGC[i]+=-para.Wgoc_distrb[i]*GoC;
	}

	// GC rates
	for (int i=0; i<para.Ngc; i++) {
		nuGC[i]=para.gain_vector[i]*(*para.Ptr_gcIO) (hGC[i]-para.threshold_vector[i],para.sp_alpha);

		// add GC interference signal
		if(para.intrfflag==2){
			for (int j=0; j<para.Ncomp_intf; j++)	nuGC[i]=nuGC[i]+input_intf[j]*para.Intf_weights.at(i,j);
		}
	}

	// GC filtering
	if(para.filterflag>0)	nuGC=filter_GCs(nuGC, para);

	//	RNN feedback input
	if(para.closedloop==0){
		// open loop: drive network with target signal
		feedback=func_target(t, para);
	}
	else {
		// closed loop: drive network with output signal
		feedback=0;
		for(int i= 0; i < para.Ngc; i++) {
			feedback+=nuGC[i]*para.Wout_PC[i];
		}
		feedback=feedback/double(para.Ngc);
	}

	// add RNN interference signal
	if(para.intrfflag==1){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
		for(int i= 0; i < para.Ncx; i++) {
			interference[i]=0;
			for (int j=0; j<para.Ncomp_intf; j++){
				interference[i]+=input_intf[j]*para.Wfb_interf.at(i,j);
			}
		}
	}

	// recurrent inputs
	for(int i= 0; i < para.Ncx; i++) {
		recurrent[i]=0;
		for (int j=0; j<para.Wcx[i].Nsyn; j++) recurrent[i]+=nuCx[para.Wcx[i].idx[j]]*para.Wcx[i].weight[j];
	}

	for(int i= 0; i < para.Ncx; i++) {
		dhcx_dt[i]= (-hcx[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	return GSL_SUCCESS;
}
