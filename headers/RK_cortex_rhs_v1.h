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
