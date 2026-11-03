int RK_Floquet_cortex_rhs(double t, const double X[], double dX_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);
	arma::vec dPHI(para.Ncx);
	arma::vec nuCx(para.Ncx);

	// extract state variables from X array
	const arma::mat Sm(const_cast<double*>(&X[para.Ncx]),para.Ncx,para.Ncx,false,true);

	// firing rates (no noise)
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	//	feedback input
	// open loop: drive network with target signal
	if(para.floquetmodeflag==0) {
		feedback=func_target(t, para);
	}
	// closed loop: drive network with output signal
	else if(para.floquetmodeflag==1) {
		feedback=0;
		for(int i= 0; i < para.Ncx; i++) {
			feedback+=nuCx[i]*para.Wout[i];
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
		for (int j=0; j<para.Wcx[i].Nsyn; j++) recurrent[i]+=nuCx[para.Wcx[i].idx[j]]*para.Wcx[i].weight[j];
	}

	for(int i= 0; i < para.Ncx; i++) {
		dX_dt[i]= (-X[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	// *************************
	// update monodromy matrix
	// *************************
	arma::mat dSm_dt;

	// tanh derivative
	if(para.actfuncflag==0) {
		for(int i= 0; i < para.Ncx; i++) {
			dPHI[i]=1.0 - nuCx[i]*nuCx[i];
		}
	}

	// linear RNN with or without interference signal
	if(para.actfuncflag==1){
		dSm_dt=para.Weff *Sm;
	}
	// tanh RNN
	else if(para.actfuncflag==0){
		arma::mat Weff=-para.I+para.WrankN*arma::diagmat(dPHI);
		dSm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *Sm;
	}

	// assign monodromy matrix to dX_dt array
	std::memcpy(dX_dt + para.Ncx, dSm_dt.memptr(), sizeof(double) * para.Ncx*para.Ncx);

	return GSL_SUCCESS;
}

int RK_Floquet_cortex_cereb_rhs(double t, const double X[], double dX_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	double GoC;
	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);
	arma::vec MF(para.Nmf);
	arma::vec dPHI(para.Ncx);
	arma::vec dPSI(para.Ngc);
	arma::vec nuCx(para.Ncx);
	arma::vec nuGC(para.Ngc);
	arma::vec hGC(para.Ngc);

	// extract state variables from X array
	const arma::mat Sm(const_cast<double*>(&X[para.Ncx]),para.Ncx,para.Ncx,false,true);

	// RNN firing rates (no noise)
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// GC interference signal
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
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

	// Instantaneous GCs
	for (int i=0; i<para.Ngc; i++) {
		//	pre-synaptic inputs of GC i
		double inputs=0;
		for (int j=0; j<para.Nsyn; j++) {
			inputs+=MF[para.Wmf[i].idx[j]];
			//~ inputs+=Wmf[i].weight[j]*MF[Wmf[i].idx[j]];
		}
		inputs=inputs/para.sqNsyn;
		if(para.gocmodeflag>0) inputs+=-para.Wgoc_distrb[i]*GoC;

		// instananeous GCs
		hGC[i]=inputs;
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

	//	feedback input
	// open loop: drive network with target signal
	if(para.floquetmodeflag==0) {
		feedback=func_target(t, para);
	}
	// closed loop: drive network with output signal
	else if(para.floquetmodeflag==1) {
		feedback=0;
		for(int i= 0; i < para.Ngc; i++) {
			feedback+=nuGC[i]*para.Wout_PC[i];
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
		for (int j=0; j<para.Wcx[i].Nsyn; j++) recurrent[i]+=nuCx[para.Wcx[i].idx[j]]*para.Wcx[i].weight[j];
	}

	for(int i= 0; i < para.Ncx; i++) {
		dX_dt[i]= (-X[i] + recurrent[i]+feedback*para.Wfb[i]+interference[i])/para.tau_distrb[i];
	}

	// *************************
	// update monodromy matrix
	// *************************
	arma::mat dSm_dt;

	// tanh derivative
	if(para.actfuncflag==0) {
		for(int i= 0; i < para.Ncx; i++) {
			dPHI[i]=1.0 - nuCx[i]*nuCx[i];
		}
	}
	// ReLu derivative (heaviside)
	if(para.actfuncflag_GC==0){
		for(int i= 0; i < para.Ngc; i++) {
			dPSI[i]=para.gain_vector[i]*heavi(hGC[i]-para.threshold_vector[i]);
		}
	}
	// softplus derivative
	else if(para.actfuncflag_GC==1){
		for(int i= 0; i < para.Ngc; i++) {
			dPSI[i]=para.gain_vector[i]*inv_logit(hGC[i]-para.threshold_vector[i], para.sp_alpha);
		}
	}

	// without Golgi cell feedback
	if(para.gocmodeflag==0){
		if(para.actfuncflag==1 && para.actfuncflag_GC==2){
			// linear RNN & linear GCs
			dSm_dt=para.Weff *Sm;
		}
		else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Weff=para.WrankN+para.Wrank1*arma::diagmat(dPSI)*para.Wctxgc;
			dSm_dt=Weff *Sm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Weff=-para.I+(para.Wcortex+para.Wrank1*arma::diagmat(dPSI)*para.Wctxgc)*arma::diagmat(dPHI);
			dSm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *Sm;
		}
	}
	// with Golgi cell feedback
	else if (para.gocmodeflag==1){
		arma::mat Igc= arma::eye( para.Ngc, para.Ngc );
		arma::mat Wtmp( para.Ngc, para.Ngc );
		if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Wtmp=arma::inv(Igc+para.Wgoc_distrb_scaled*arma::trans(dPSI));
			arma::mat Weff=para.WrankN+para.Wrank1*arma::diagmat(dPSI)*Wtmp*para.Wctxgc;
			dSm_dt=Weff *Sm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Wtmp=arma::inv(Igc+para.Wgoc_distrb_scaled*arma::trans(dPSI));
			arma::mat Weff=-para.I+(para.Wcortex+para.Wrank1*arma::diagmat(dPSI)*Wtmp*para.Wctxgc)*arma::diagmat(dPHI);
			dSm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *Sm;
		}
	}
	// with feedforward Golgi cell inhibition
	else if (para.gocmodeflag==3)	{
		double SD_rate_mfs=arma::stddev(nuCx);
		double avrg_rate_mfs=arma::mean(nuCx);
		if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Weff=para.WrankN +para.Wrank1*arma::diagmat(dPSI)*(para.Wctxgc -para.Wgoc_distrb_scaled*arma::trans(nuCx-avrg_rate_mfs)/SD_rate_mfs);
			dSm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *Sm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			//~ arma::mat Weff=-para.I+para.Wcortex*arma::diagmat(dPHI)+para.Wrank1*arma::diagmat(dPSI)*(para.Wctxgc*arma::diagmat(dPHI) -para.Wgoc_distrb_scaled*arma::trans( dPHI%(nuCx-avrg_rate_mfs))/SD_rate_mfs);

			arma::vec Y(para.Ncx);
			for(int i= 0; i < para.Ncx; i++) {
				Y[i]= dPHI[i]*(nuCx[i]-avrg_rate_mfs);
			}
			arma::mat Weff=-para.I+para.Wcortex*arma::diagmat(dPHI)+para.Wrank1*arma::diagmat(dPSI)*(para.Wctxgc*arma::diagmat(dPHI) -para.Wgoc_distrb_scaled*arma::trans( Y)/SD_rate_mfs);

			  // GoC Jacobian
			// dGoC/dx_j = dPHI[j] * (MF_j - mean(MF)) / ((Nmf-1) * GoC)  for j in MF indices
			//~ arma::vec dGoC_dX = arma::zeros(para.Ncx);
			//~ for (int j = 0; j < para.Nmf; j++) {
				//~ int idx = para.Ncx_to_MF_idx[j];
				//~ dGoC_dX[idx] = dPHI[idx] * (MF[j] - avrg_rate_mfs) / ((para.Nmf - 1) *SD_rate_mfs);
			//~ }
			//~ arma::mat Wtmp = para.Wgoc_distrb_scaled * arma::trans(dGoC_dX);
			//~ arma::mat Weff=-para.I+para.Wcortex*arma::diagmat(dPHI)+para.Wrank1*arma::diagmat(dPSI)*(para.Wctxgc*arma::diagmat(dPHI) -Wtmp);

			dSm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *Sm;
		}
	}

	// assign monodromy matrix to dX_dt array
	std::memcpy(dX_dt + para.Ncx, dSm_dt.memptr(), sizeof(double) * para.Ncx*para.Ncx);

	return GSL_SUCCESS;
}

int RK_Floquet_cortex_cereb_dyn_rhs(double t, const double X[], double dX_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	const int Nall=para.Ncx+para.Ngc;
	const int Nall2=Nall*Nall;
	double GoC;
	double feedback;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);
	arma::vec GCinputs(para.Ngc);
	arma::vec MF(para.Nmf);
	arma::vec dPHI(para.Ncx);
	arma::vec dPSI(para.Ngc);
	arma::vec nuCx(para.Ncx);
	arma::vec nuGC(para.Ngc);

    const arma::mat Sm(const_cast<double*>(X + Nall), Nall, Nall, false, true);

	// RNN firing rates (no noise)
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

	// GC filtering
	if(para.filterflag>0)	nuGC=filter_GCs(nuGC, para);

	//	feedback input
	// open loop: drive network with target signal
	if(para.floquetmodeflag==0) {
		feedback=func_target(t, para);
	}
	// closed loop: drive network with output signal
	else if(para.floquetmodeflag==1) {
		feedback=0;
		for(int i= 0; i < para.Ngc; i++) {
			feedback+=nuGC[i]*para.Wout_PC[i];
		}
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

	// GC inputs
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
		dX_dt[i+para.Ncx]= (-X[i+para.Ncx] + GCinputs[i])/para.tauGC_distrb[i];
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

	// *************************
	// update monodromy matrix
	// *************************
	arma::mat dSm_dt;

	// tanh derivative
	if(para.actfuncflag==0) {
		for(int i= 0; i < para.Ncx; i++) {
			dPHI[i]=1.0 - nuCx[i]*nuCx[i];
		}
	}

	// ReLu derivative (heaviside)
	if(para.actfuncflag_GC==0){
		for(int i= 0; i < para.Ngc; i++) {
			dPSI[i]=para.gain_vector[i]*heavi(X[i+para.Ncx]-para.threshold_vector[i]);
		}
	}
	// softplus derivative
	else if(para.actfuncflag_GC==1){
		for(int i= 0; i < para.Ngc; i++) {
			dPSI[i]=para.gain_vector[i]*inv_logit(X[i+para.Ncx]-para.threshold_vector[i], para.sp_alpha);
		}
	}

	// without Golgi cell feedback
	if(para.gocmodeflag==0){
		if(para.actfuncflag==1 && para.actfuncflag_GC==2){
			// linear RNN & cereb linear
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*para.Weff *Sm;
		}
		else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1*arma::diagmat(dPSI));
			arma::mat Weff=arma::join_vert( AB, para.CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat AB = arma::join_horiz( -para.I + para.Wcortex*arma::diagmat(dPHI), para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc*arma::diagmat(dPHI),-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
	}
	// with Golgi cell feedback
	else if (para.gocmodeflag==1){
		if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Igc+para.Wgoc_distrb_scaled*arma::trans(dPSI);
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc,-Wtmp);
			arma::mat Weff=arma::join_vert( AB, CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Igc+para.Wgoc_distrb_scaled*arma::trans(dPSI);
			arma::mat AB = arma::join_horiz( -para.I+para.Wcortex*arma::diagmat(dPHI), para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc*arma::diagmat(dPHI),-Wtmp);
			arma::mat Weff=arma::join_vert( AB, CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
	}
	else if (para.gocmodeflag==3){
		double SD_rate_mfs=arma::stddev(MF);
		double avrg_rate_mfs=arma::mean(MF);
		if(para.actfuncflag==1 && para.actfuncflag_GC==2){
			// linear RNN & cereb linear
			arma::mat Wtmp = para.Wgoc_distrb_scaled*arma::trans(nuCx-avrg_rate_mfs)/SD_rate_mfs;
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1);
			arma::mat CD = arma::join_horiz(para.Wctxgc-Wtmp,-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
		else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Wgoc_distrb_scaled*arma::trans(nuCx-avrg_rate_mfs)/SD_rate_mfs;
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc-Wtmp,-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Wgoc_distrb_scaled*arma::trans( dPHI%(nuCx-avrg_rate_mfs))/SD_rate_mfs;
			arma::mat AB = arma::join_horiz( -para.I+para.Wcortex*arma::diagmat(dPHI), para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc*arma::diagmat(dPHI)-Wtmp,-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dSm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *Sm;
		}
	}

	// assign monodromy matrix to dX_dt array
	std::memcpy(dX_dt + Nall, dSm_dt.memptr(), sizeof(double) * Nall2);


	return GSL_SUCCESS;
}

void calc_floquet_RK(double XIN[], netpara &para, nlohmann::json &json_measures, std::ofstream &alg_steps_file) {

if(para.floquetflag==1){
	if(para.algswitch==1)	alg_steps_file << "\t \t \t Calculate Floquet multipliers for ctx"<< std::endl;
	else					alg_steps_file << "\t \t \t Calculate Floquet multipliers for ctx + cereb"<< std::endl;

	para.inv_tau_distrb=1.0/para.tau_distrb;
	para.inv_tauGC_distrb=1.0/para.tauGC_distrb;

	// recurrent weights
	para.Wcortex= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
	for(int i= 0; i < para.Ncx; i++) {
		for (int j=0; j<para.Wcx[i].Nsyn; j++) {
			para.Wcortex.at(i,para.Wcx[i].idx[j])=para.Wcx[i].weight[j];
		}
	}
	arma::mat Sm;
	para.I = arma::eye( para.Ncx, para.Ncx );

	// 	set up GSL integrator
	size_t nn2;
	size_t nnet;
	gsl_odeiv2_system FLOQsys;
	gsl_odeiv2_system NETsys;

	// *************************************
	// set up matrices for cortex only
	// *************************************
	//~ if(para.cortexonlyflag==1 && para.algswitch==1){		// only compute FM when in CORTEXonly mode
	if(para.algswitch==1){
		Sm = arma::eye( para.Ncx, para.Ncx );
		para.Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
		// linear RNN with interference signal
		if(para.actfuncflag==1){
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					para.Weff.at(i,j)=-para.I.at(i,j)+para.Wcortex.at(i,j) + para.Wfb[i]*para.Wout[j];
				}
			}
			para.Weff= arma::diagmat(para.inv_tau_distrb)*para.Weff;
		}
		// tanh RNN
		else if(para.actfuncflag==0){
			para.WrankN= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					para.WrankN.at(i,j)=para.Wcortex.at(i,j) + para.Wfb[i]*para.Wout[j];
				}
			}
		}
		nn2 = static_cast<size_t>(para.Ncx+para.Ncx*para.Ncx);
		FLOQsys = {RK_Floquet_cortex_rhs, NULL, nn2, &para};

		nnet = static_cast<size_t>(para.Ncx);
		NETsys = {RK_cortex_rhs, NULL, nnet, &para};

	}

	// *************************************
	// set up matrices for cortex + cerebellum
	// *************************************
	if(para.algswitch==2){
		// MF-GC weights
		arma::mat Wmfgc= arma::zeros<arma::mat>(para.Ngc,para.Nmf);
		for(int i= 0; i < para.Ngc; i++) {
			for (int j=0; j<para.Nsyn; j++) {
				Wmfgc.at(i,para.Wmf[i].idx[j])=1.0/para.sqNsyn;
				//~ Wmfgc.at(i,para.Wmf[i].idx[j])=1.0;
			}
		}
		// CTX-MF weights
		arma::mat Wctxmf= arma::zeros<arma::mat>(para.Nmf,para.Ncx);
		for(int i= 0; i < para.Nmf; i++) {
			//~ Wctxmf.at(i,para.Ncx_to_MF_idx[i])=1.0;
			Wctxmf.at(i,para.Ncx_to_MF_idx[i])=1.0;
		}
		// effective CTX-GC weights
		para.Wctxgc=Wmfgc*Wctxmf;

		// rank 1 matrix
		para.Wrank1= arma::zeros<arma::mat>(para.Ncx,para.Ngc);
		for(int i= 0; i < para.Ncx; i++) {
			for (int j=0; j<para.Ngc; j++) {
				para.Wrank1.at(i,j)= para.Wfb[i]*para.Wout_PC[j];
			}
		}

		// scaled GoC weights
		para.Wgoc_distrb_scaled=para.Wgoc_distrb/double(para.Ngc);

		// *************************************
		// instantaneous GCs
		// *************************************
		if(para.taugcmodeflag==0) {
			Sm = arma::eye( para.Ncx, para.Ncx );
			para.Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);

			// rank N matrix
			para.WrankN=(-para.I+para.Wcortex)*arma::diagmat(para.inv_tau_distrb);

			// without Golgi cell feedback
			if(para.gocmodeflag==0){
				if(para.actfuncflag==1 && para.actfuncflag_GC==2){
					// linear RNN & cereb
					para.Wrank1=para.Wrank1*arma::diagmat(para.gain_vector)*para.Wctxgc*arma::diagmat(para.inv_tau_distrb);
					para.Weff=para.WrankN+para.Wrank1;
				}
				else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					para.Wctxgc=para.Wctxgc*arma::diagmat(para.inv_tau_distrb);
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
			// with Golgi cell feedback
			else if (para.gocmodeflag==1){
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					para.Wctxgc=para.Wctxgc*arma::diagmat(para.inv_tau_distrb);
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
			else if (para.gocmodeflag==2)	std::cerr << "Heterogeneous Golgi cell feedback not implemented (yet)." << std::endl;
			// with feedforward Golgi cell inhibition
			else if (para.gocmodeflag==3)	{
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					para.WrankN=(-para.I+para.Wcortex);
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}

			nn2 = static_cast<size_t>(para.Ncx+para.Ncx*para.Ncx);
			FLOQsys = {RK_Floquet_cortex_cereb_rhs, NULL, nn2, &para};

			nnet = static_cast<size_t>(para.Ncx);
			NETsys = {RK_cortex_cereb_rhs, NULL, nnet, &para};
		}
		// *************************************
		// non-instantaneous GCs
		// *************************************
		else{
			Sm = arma::eye( para.Ncx + para.Ngc, para.Ncx + para.Ngc);
			para.Weff=arma::zeros<arma::mat>( para.Ncx + para.Ngc, para.Ncx + para.Ngc);
			para.inv_tauALL_distrb = arma::join_vert(para.inv_tau_distrb,para.inv_tauGC_distrb);
			// rank N matrix
			para.WrankN=-para.I + para.Wcortex;
			para.Igc= arma::eye( para.Ngc, para.Ngc );

			// without Golgi cell feedback
			if(para.gocmodeflag==0){
				if(para.actfuncflag==1 && para.actfuncflag_GC==2){
					// linear RNN & cereb
					arma::mat AB = arma::join_horiz(para.WrankN,para.Wrank1);
					arma::mat CD = arma::join_horiz(para.Wctxgc,-para.Igc);
					para.Weff = arma::join_vert(AB,CD);
				}
				else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					para.CD = arma::join_horiz(para.Wctxgc,-para.Igc);
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
			// with Golgi cell feedback
			else if (para.gocmodeflag==1){
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;

			}
			else if (para.gocmodeflag==2)	std::cerr << "Heterogeneous Golgi cell feedback not implemented (yet)." << std::endl;
			else if (para.gocmodeflag==3){
				if(para.actfuncflag==1 && para.actfuncflag_GC==2){
					// linear RNN & cereb
				}
				else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}

			nn2 = static_cast<size_t>(para.Ncx+para.Ngc+(para.Ncx+para.Ngc)*(para.Ncx+para.Ngc));
			FLOQsys = {RK_Floquet_cortex_cereb_dyn_rhs, NULL, nn2, &para};

			nnet = static_cast<size_t>(para.Ncx+para.Ngc);
			NETsys = {RK_cortex_cereb_dyn_rhs, NULL, nnet, &para};
		}
	}

	// *****************************
	// integrate system
	// *****************************
	//~ gsl_odeiv2_driver * FLOQdriver = gsl_odeiv2_driver_alloc_y_new (&FLOQsys, gsl_odeiv2_step_rk8pd,para.dt, para.INT_delta_abs,para.INT_delta_rel);
	gsl_odeiv2_driver * FLOQdriver = gsl_odeiv2_driver_alloc_y_new (&FLOQsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs/10.0,para.INT_delta_rel/10.0);

	gsl_odeiv2_driver * NETdriver = gsl_odeiv2_driver_alloc_y_new (&NETsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs/10.0,para.INT_delta_rel/10.0);

	std::vector<double> Xnet(nnet);
	std::vector<double> X(nn2);
	int status;
	double time=0;

	// *****************************
	// Closed-loop transient period
	// *****************************
	// assign state variables
	if(para.taugcmodeflag==0){
		for(int i= 0; i < para.Ncx; i++) {
			Xnet[i]=XIN[i];
		}
	}
	else if(para.taugcmodeflag>0){
		for(int i= 0; i < para.Ncx+para.Ngc; i++) {
			Xnet[i]=XIN[i];
		}
	}

	para.closedloop=1;		// closed loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t FM transient phase" << "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	status = gsl_odeiv2_driver_apply (NETdriver, &time, (para.Tmax_trans-1)*para.dt, Xnet.data());
	if (status != GSL_SUCCESS){
		printf ("Error during FM transient phase, return value=%d\n", status);
	}
	alg_steps_file << "\t \t \t End of FM transient phase, T =  "<< para.Tmax_trans*para.dt << "\t Ttot =  "<< time << std::endl;

	// *****************************
	// Integrate monodromy matrix for one period
	// *****************************
	// assign state variables
	if(para.taugcmodeflag==0){
		for(int i= 0; i < para.Ncx; i++) {
			X[i]=Xnet[i];
			for(int j= 0; j < para.Ncx; j++) {
				int k=para.Ncx+i+j*para.Ncx;
				X[k]=Sm.at(i,j);
			}
		}
	}
	else if(para.taugcmodeflag>0){
		for(int i= 0; i < para.Ncx+para.Ngc; i++) {
			X[i]=Xnet[i];
		}
		for(int i= 0; i < para.Ncx+para.Ngc; i++) {
			for(int j= 0; j < para.Ncx+para.Ngc; j++) {
				int k=para.Ncx+para.Ngc+i+j*(para.Ncx+para.Ngc);
				X[k]=Sm.at(i,j);
			}
		}
	}
	//~ std::cout << arma::size(Sm) << std::endl;

	time=0;
	para.closedloop=1;		// closed loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t FM computation phase" << "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	// full monodromy matrix
	status = gsl_odeiv2_driver_apply (FLOQdriver, &time, (para.Tmax_FM)*para.dt, X.data());
	if (status != GSL_SUCCESS){
		printf ("Error during FM phase, return value=%d\n", status);
	}
	alg_steps_file << "\t \t \t End of FM computation phase, T =  "<< para.Tmax_FM*para.dt << "\t Ttot =  "<< time << std::endl;

	// extract state variables from X array
	if(para.taugcmodeflag==0){
		for(int i= 0; i < para.Ncx; i++) {
			for(int j= 0; j < para.Ncx; j++) {
				int k=para.Ncx+i+j*para.Ncx;
				Sm.at(i,j)=X[k];
			}
		}
	}
	else if(para.taugcmodeflag>0){
		for(int i= 0; i < para.Ncx+para.Ngc; i++) {
			for(int j= 0; j < para.Ncx+para.Ngc; j++) {
				int k=para.Ncx+para.Ngc+i+j*(para.Ncx+para.Ngc);
				Sm.at(i,j)=X[k];
			}
		}
	}

	// calculate EVs of monodromy matrix Sm
	int NEV=para.Ncx;
	if(para.taugcmodeflag>0) NEV=para.Ncx+para.Ngc;
	arma::cx_vec Sm_eigval(NEV);
	if(Sm.has_nan()==0){
		Sm_eigval = arma::eig_gen( Sm );
		Sm_eigval=Sm_eigval(arma::sort_index( arma::abs(Sm_eigval), "descend" ));
	}
	else 	Sm_eigval(0)=arma::datum::nan;
	//~ std::cout << arma::size(Sm_eigval) << std::endl;

	if(para.learnflag<99){
			if(para.algswitch==1)	std::cout << "RNN Floquet multipliers:" << std::endl;
			else 					std::cout << "Cereb Floquet multipliers:" << std::endl;
			// for(int i= 0; i < 4; i++) std::cout << std::real(Sm_eigval[i]) <<  "\t" << std::imag(Sm_eigval[i]) <<  "\t" << std::abs(Sm_eigval[i]) << std::endl;
			for(int i= 0; i < 4; i++) std::cout << std::abs(Sm_eigval[i]) << std::endl;
			std::cout << std::endl;
	}

	for(int i= 0; i < 1; i++) {
		if(para.algswitch==1)	json_measures["Max_Floquet_multiplier_ctx"].push_back(std::abs(Sm_eigval[i]));
		else					json_measures["Max_Floquet_multiplier_cereb"].push_back(std::abs(Sm_eigval[i]));
	}

	gsl_odeiv2_driver_free (FLOQdriver);
	gsl_odeiv2_driver_free (NETdriver);

}
else{
	if(para.algswitch==1)	json_measures["Max_Floquet_multiplier_ctx"].push_back(NAN);
	else					json_measures["Max_Floquet_multiplier_cereb"].push_back(NAN);
}

	return;
}
