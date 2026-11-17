////////////////////////////////////////////////////
// RNN only
////////////////////////////////////////////////////
int RK_cortex_rhs(double t, const double hcx[], double dhcx_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);

	// firing rates
	for(int i= 0; i < para.Ncx; i++) {
		para.nuCx[i]=(*para.Ptr_ctxIO)(hcx[i],1.0);
	}

	//	feedback input
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
		feedback=feedback/double(para.Ncx);
	}

	// add interference signal
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
		dhcx_dt[i]= (-hcx[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	return GSL_SUCCESS;
}

////////////////////////////////////////////////////
// RNN + instantaneous FFL
////////////////////////////////////////////////////
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
	if(para.gocmodeflag==1){
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

////////////////////////////////////////////////////
// RNN + dynamic FFL
////////////////////////////////////////////////////
int RK_cortex_cereb_dyn_rhs(double t, const double X[], double dX_dt[], void *params) {

	// hcx are the first Nctx entries of X
	// hgc are the last Ngc entries of X

	auto& para = *static_cast<netpara*>(params);

	double GoC;
	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference= arma::zeros<arma::vec>(para.Ncx);
	arma::vec hGC(para.Ngc);
	arma::vec MF(para.Nmf);
	arma::vec nuCx(para.Ncx);
	arma::vec nuGC(para.Ngc);

	// RNN firing rates (no noise)
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// MFs = sample from cortical activity
	for (int j=0; j<para.Nmf; j++) 	MF[j]=nuCx[para.Ncx_to_MF_idx[j]];

	// Golgi cell
	if(para.gocmodeflag==1){
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

	// GC inputs
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
	for(int i= 0; i < para.Ngc; i++) {
		dX_dt[i+para.Ncx]= (-X[i+para.Ncx] + hGC[i])/para.tauGC_distrb[i];
	}

	// GC rates
	for (int i=0; i<para.Ngc; i++) {
		nuGC[i]=para.gain_vector[i]*(*para.Ptr_gcIO) (X[i+para.Ncx]-para.threshold_vector[i],para.sp_alpha);

		// add GC interference signal
		if(para.intrfflag==2){
			for (int j=0; j<para.Ncomp_intf; j++)	nuGC[i]=nuGC[i]+input_intf[j]*para.Intf_weights.at(i,j);
		}
	}

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
		dX_dt[i]= (-X[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	return GSL_SUCCESS;
}
