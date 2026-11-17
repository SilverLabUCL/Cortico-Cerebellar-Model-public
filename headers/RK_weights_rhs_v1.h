////////////////////////////////////////////////////
// RNN only + DELTA weight update
////////////////////////////////////////////////////
int RK_cortex_weights_rhs(double t, const double X[], double dX_dt[], void *params) {

	// DELTA learning: update both hcx and Wout
	// hcx are the first Nctx entries of X
	// Wout are the last Nctx entries of X

	auto& para = *static_cast<netpara*>(params);

	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference= arma::zeros<arma::vec>(para.Ncx);
	arma::vec nuCx(para.Ncx);

	// firing rates
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// feedback input
	// closed loop: drive network with output signal
	feedback=0;
	for(int i= 0; i < para.Ncx; i++) {
		feedback+=nuCx[i]*X[i+para.Ncx];
	}
	feedback=feedback/double(para.Ncx);

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
		for (int j=0; j<para.Wcx[i].Nsyn; j++) recurrent[i]+=nuCx[para.Wcx[i].idx[j]]*para.Wcx[i].weight[j];
	}

	for(int i= 0; i < para.Ncx; i++) {
		dX_dt[i]= (-X[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	// update the error of the linear readout unit
	para.error_signed_train = feedback-func_target(t, para);

	// norm of firing rates
	double r2 = para.delta_r2_reg;
	for (int i=0; i<para.Ncx; i++)	r2 += nuCx[i] * nuCx[i];
	r2=r2/double(para.Ncx);

	// update the output weights
	//~ double lr_normalized = para.delta_eta / (para.tau_avrg * r2);
	double lr_normalized = para.delta_eta *para.lambda/ r2;
	double eta_error = lr_normalized * para.error_signed_train;

	for (int i=0; i<para.Ncx; i++){
		dX_dt[para.Ncx+i]=-nuCx[i]*eta_error;
	}

	return GSL_SUCCESS;
}

////////////////////////////////////////////////////
// RNN + instantaneous GCs + DELTA weight update
////////////////////////////////////////////////////
int RK_cortex_cereb_weights_rhs(double t, const double hcxW[], double dhcxW_dt[], void *params) {

	// DELTA learning: update both hcx and Wout
	// hcx are the first Nctx entries of hcxW
	// Wout are the last Nctx entries of hcxW

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
		nuCx[i]=(*para.Ptr_ctxIO)(hcxW[i],1.0);
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

	// feedback input
	// closed loop: drive network with output signal
	feedback=0;
	for(int i= 0; i < para.Ngc; i++) {
		feedback+=nuGC[i]*hcxW[i+para.Ncx];
	}
	feedback=feedback/double(para.Ngc);

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
		dhcxW_dt[i]= (-hcxW[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	// update the error of the linear readout unit
	para.error_signed_train = feedback-func_target(t, para);

	// norm of firing rates
	double r2 = para.delta_r2_reg;
	for (int i=0; i<para.Ngc; i++)	r2 += nuGC[i] * nuGC[i];
	r2=r2/double(para.Ngc);

	// update the output weights
	double lr_normalized = para.delta_eta *para.lambda/ r2;
	double eta_error = lr_normalized * para.error_signed_train;

	for (int i=0; i<para.Ngc; i++){
		dhcxW_dt[para.Ncx+i]=-nuGC[i]*eta_error;
	}

	return GSL_SUCCESS;
}

////////////////////////////////////////////////////
// RNN + dynamic GCs + DELTA weight update
////////////////////////////////////////////////////
int RK_cortex_cereb_dyn_weights_rhs(double t, const double X[], double dX_dt[], void *params) {

	// DELTA learning: update both hcx and Wout
	// hcx are the first Nctx entries of X
	// hgc are the next Ngc entries of X
	// Wout are the last Ngc entries of X

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

	// RNN firing rates
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// GC interference signal
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
	}

	// GC rates
	for (int i=0; i<para.Ngc; i++) {
		nuGC[i]=para.gain_vector[i]*(*para.Ptr_gcIO) (X[i+para.Ncx]-para.threshold_vector[i],para.sp_alpha);

		// add GC interference signal
		if(para.intrfflag==2){
			for (int j=0; j<para.Ncomp_intf; j++)	nuGC[i]=nuGC[i]+input_intf[j]*para.Intf_weights.at(i,j);
		}
	}

	// feedback input
	// closed loop: drive network with output signal
	feedback=0;
	for(int i= 0; i < para.Ngc; i++) {
		feedback+=nuGC[i]*X[i+para.Ncx+para.Ngc];
	}
	feedback=feedback/double(para.Ngc);

	// MFs = sample from cortical activity
	for (int j=0; j<para.Nmf; j++) 	MF[j]=nuCx[para.Ncx_to_MF_idx[j]];

	// Golgi cell
	if(para.gocmodeflag==1){
		GoC=arma::mean(nuGC);
	}
	else if(para.gocmodeflag==3){
		GoC=arma::stddev(MF);
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

	// update the error of the linear readout unit
	para.error_signed_train = feedback-func_target(t, para);

	// norm of firing rates
	double r2 = para.delta_r2_reg;
	for (int i=0; i<para.Ngc; i++)	r2 += nuGC[i] * nuGC[i];
	r2=r2/double(para.Ngc);

	// update the output weights
	double lr_normalized = para.delta_eta *para.lambda/ r2;
	double eta_error = lr_normalized * para.error_signed_train;

	for (int i=0; i<para.Ngc; i++){
		dX_dt[para.Ncx+para.Ngc+i]=-nuGC[i]*eta_error;
	}

	return GSL_SUCCESS;
}

