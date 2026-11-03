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
	double r2 = 1e-6;
	for (int i=0; i<para.Ncx; i++)	r2 += nuCx[i] * nuCx[i];
	r2=r2/double(para.Ncx);

	// update the output weights
	//~ double lr_normalized = para.delta_eta / (para.tau_avrg * r2);
	double lr_normalized = para.delta_eta *para.lambda*10/ r2;
	double eta_error = lr_normalized * para.error_signed_train;

	for (int i=0; i<para.Ncx; i++){
		dX_dt[para.Ncx+i]=-nuCx[i]*eta_error;//  -para.delta_eps_decay * X[para.Ncx+i];
	}

	return GSL_SUCCESS;
}
