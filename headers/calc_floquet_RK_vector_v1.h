int RK_Floquet_cortex_rhs_vector(double t, const double X[], double dX_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	double feedback;
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);
	arma::vec dPHI(para.Ncx);
	arma::vec nuCx(para.Ncx);

	// extract state variables from X array
	const arma::vec vm(X + para.Ncx, para.Ncx);

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
		feedback=feedback/double(para.Ncx);
	}

	// add interference signal
	if(para.intrfflag==1){
		arma::vec input_intf(para.Ncomp_intf);
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
	arma::vec dvm_dt;

	// tanh derivative
	if(para.actfuncflag==0) {
		for(int i= 0; i < para.Ncx; i++) {
			dPHI[i]=1.0 - nuCx[i]*nuCx[i];
		}
	}

	// linear RNN with or without interference signal
	if(para.actfuncflag==1){
		dvm_dt=para.Weff *vm;
	}
	// tanh RNN
	else if(para.actfuncflag==0){
		arma::mat Weff=-para.I+para.WrankN*arma::diagmat(dPHI);
		dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
	}

	// assign vm vector to dX_dt array
	std::memcpy(dX_dt + para.Ncx, dvm_dt.memptr(), sizeof(double) * para.Ncx);

	return GSL_SUCCESS;
}

int RK_Floquet_cortex_cereb_rhs_vector(double t, const double X[], double dX_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	double GoC;
	double feedback;
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);
	arma::vec MF(para.Nmf);
	arma::vec dPHI(para.Ncx);
	arma::vec dPSI(para.Ngc);
	arma::vec nuCx(para.Ncx);
	arma::vec nuGC(para.Ngc);
	arma::vec hGC(para.Ngc);

	// extract state variables from X array
	const arma::vec vm(X + para.Ncx, para.Ncx);

	// RNN firing rates (no noise)
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// GC interference signal
	arma::vec input_intf(para.Ncomp_intf);
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
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

	// *************************
	// update monodromy matrix
	// *************************
	arma::vec dvm_dt;

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
			dvm_dt=arma::diagmat(para.inv_tau_distrb)*para.Weff *vm;
		}
		else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Weff=para.WrankN+para.Wrank1*arma::diagmat(dPSI)*para.Wctxgc;
			dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Weff=-para.I+(para.Wcortex+para.Wrank1*arma::diagmat(dPSI)*para.Wctxgc)*arma::diagmat(dPHI);
			dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
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
			dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Wtmp=arma::inv(Igc+para.Wgoc_distrb_scaled*arma::trans(dPSI));
			arma::mat Weff=-para.I+(para.Wcortex+para.Wrank1*arma::diagmat(dPSI)*Wtmp*para.Wctxgc)*arma::diagmat(dPHI);
			dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
		}
	}
	// with feedforward Golgi cell inhibition (tested for Nmf=Ncx)
	else if (para.gocmodeflag==3)	{
		double SD_rate_mfs=GoC;
		double avrg_rate_mfs=arma::mean(MF);
		double SD_rate_mfs_safe = std::max(GoC, 1e-10);
		if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Weff=para.WrankN +para.Wrank1*arma::diagmat(dPSI)*(para.Wctxgc -para.Wgoc_distrb_scaled*arma::trans(nuCx-avrg_rate_mfs)/SD_rate_mfs_safe);
			dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::vec Y(para.Ncx);
			for(int i= 0; i < para.Ncx; i++) {
				Y[i]= dPHI[i]*(nuCx[i]-avrg_rate_mfs);
			}
			arma::mat Wtmp = para.Wgoc_distrb_scaled * arma::trans(Y) / SD_rate_mfs_safe;
			arma::mat Weff=-para.I+para.Wcortex*arma::diagmat(dPHI)+para.Wrank1*arma::diagmat(dPSI)*(para.Wctxgc*arma::diagmat(dPHI) -Wtmp);

			dvm_dt=arma::diagmat(para.inv_tau_distrb)*Weff *vm;
		}
	}

	// assign vm vector to dX_dt array
	std::memcpy(dX_dt + para.Ncx, dvm_dt.memptr(), sizeof(double) * para.Ncx);

	return GSL_SUCCESS;
}

int RK_Floquet_cortex_cereb_dyn_rhs_vector(double t, const double X[], double dX_dt[], void *params) {

	auto& para = *static_cast<netpara*>(params);

	const int Nall=para.Ncx+para.Ngc;
	double GoC;
	double feedback;
	arma::vec recurrent(para.Ncx);
	arma::vec interference(para.Ncx, arma::fill::zeros);
	arma::vec GCinputs(para.Ngc);
	arma::vec MF(para.Nmf);
	arma::vec dPHI(para.Ncx);
	arma::vec dPSI(para.Ngc);
	arma::vec nuCx(para.Ncx);
	arma::vec nuGC(para.Ngc);

	// extract state variables from X array
	const arma::vec vm(X + Nall, Nall);

	// RNN firing rates (no noise)
	for(int i= 0; i < para.Ncx; i++) {
		nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// GC interference signal
	arma::vec input_intf(para.Ncomp_intf);
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
		feedback=feedback/double(para.Ngc);
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
	arma::vec dvm_dt;

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
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*para.Weff *vm;
		}
		else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1*arma::diagmat(dPSI));
			arma::mat Weff=arma::join_vert( AB, para.CD);
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat AB = arma::join_horiz( -para.I + para.Wcortex*arma::diagmat(dPHI), para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc*arma::diagmat(dPHI),-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
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
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Igc+para.Wgoc_distrb_scaled*arma::trans(dPSI);
			arma::mat AB = arma::join_horiz( -para.I+para.Wcortex*arma::diagmat(dPHI), para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc*arma::diagmat(dPHI),-Wtmp);
			arma::mat Weff=arma::join_vert( AB, CD);
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
		}
	}
	else if (para.gocmodeflag==3){
		double SD_rate_mfs=GoC;
		double avrg_rate_mfs=arma::mean(MF);
		double SD_rate_mfs_safe = std::max(GoC, 1e-10);
		if(para.actfuncflag==1 && para.actfuncflag_GC==2){
			// linear RNN & cereb linear
			arma::mat Wtmp = para.Wgoc_distrb_scaled*arma::trans(nuCx-avrg_rate_mfs)/SD_rate_mfs_safe;
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1);
			arma::mat CD = arma::join_horiz(para.Wctxgc-Wtmp,-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
			// linear RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Wgoc_distrb_scaled*arma::trans(nuCx-avrg_rate_mfs)/SD_rate_mfs_safe;
			arma::mat AB = arma::join_horiz( para.WrankN, para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc-Wtmp,-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
		}
		else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
			// tanh RNN & cereb ReLu or softplus
			arma::mat Wtmp = para.Wgoc_distrb_scaled*arma::trans( dPHI%(nuCx-avrg_rate_mfs))/SD_rate_mfs_safe;
			arma::mat AB = arma::join_horiz( -para.I+para.Wcortex*arma::diagmat(dPHI), para.Wrank1*arma::diagmat(dPSI));
			arma::mat CD = arma::join_horiz(para.Wctxgc*arma::diagmat(dPHI)-Wtmp,-para.Igc);
			arma::mat Weff=arma::join_vert( AB, CD);
			dvm_dt=arma::diagmat(para.inv_tauALL_distrb)*Weff *vm;
		}
	}

	// assign vm vector to dX_dt array
	std::memcpy(dX_dt + Nall, dvm_dt.memptr(), sizeof(double) * Nall);

	return GSL_SUCCESS;
}

void calc_floquet_RK_vector(double XIN[], netpara &para, nlohmann::json &json_measures, std::ofstream &log_file) {

if(para.floquetflag==1){
	if(para.algswitch==1)	log_file << "\t \t \t Calculate Floquet multipliers for ctx"<< std::endl;
	else					log_file << "\t \t \t Calculate Floquet multipliers for ctx + cereb"<< std::endl;

	para.inv_tau_distrb=1.0/para.tau_distrb;
	para.inv_tauGC_distrb=1.0/para.tauGC_distrb;

	// recurrent weights
	para.Wcortex= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
	for(int i= 0; i < para.Ncx; i++) {
		for (int j=0; j<para.Wcx[i].Nsyn; j++) {
			para.Wcortex.at(i,para.Wcx[i].idx[j])=para.Wcx[i].weight[j];
		}
	}
	arma::vec vm;
	para.I = arma::eye( para.Ncx, para.Ncx );

	// 	set up GSL integrator
	size_t nn2;
	size_t nnet;
	gsl_odeiv2_system FLOQsys;

	// *************************************
	// set up matrices for cortex only
	// *************************************
	if(para.algswitch==1){

		para.Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
		// linear RNN with interference signal
		if(para.actfuncflag==1){
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					para.Weff.at(i,j)=-para.I.at(i,j)+para.Wcortex.at(i,j) + para.Wfb[i]*para.Wout[j]/double(para.Ncx);
				}
			}
			para.Weff= arma::diagmat(para.inv_tau_distrb)*para.Weff;
		}
		// tanh RNN
		else if(para.actfuncflag==0){
			para.WrankN= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					para.WrankN.at(i,j)=para.Wcortex.at(i,j) + para.Wfb[i]*para.Wout[j]/double(para.Ncx);
				}
			}
		}
		nnet = static_cast<size_t>(para.Ncx);
		nn2 = static_cast<size_t>(para.Ncx+para.Ncx);
		FLOQsys = {RK_Floquet_cortex_rhs_vector, NULL, nn2, &para};
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
			}
		}
		// CTX-MF weights
		arma::mat Wctxmf= arma::zeros<arma::mat>(para.Nmf,para.Ncx);
		for(int i= 0; i < para.Nmf; i++) {
			Wctxmf.at(i,para.Ncx_to_MF_idx[i])=1.0;
		}
		// effective CTX-GC weights
		para.Wctxgc=Wmfgc*Wctxmf;

		// rank 1 matrix
		para.Wrank1= arma::zeros<arma::mat>(para.Ncx,para.Ngc);
		for(int i= 0; i < para.Ncx; i++) {
			for (int j=0; j<para.Ngc; j++) {
				para.Wrank1.at(i,j)= para.Wfb[i]*para.Wout_PC[j]/double(para.Ngc);
			}
		}

		// scaled GoC weights
		para.Wgoc_distrb_scaled=para.Wgoc_distrb/double(para.Ncx); // normalised by Ncx; this term comes from the derivative w/r to the RNN inputs

		// *************************************
		// instantaneous GCs
		// *************************************
		if(para.taugcmodeflag==0) {
			para.Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);

			// rank N matrix
			para.WrankN=(-para.I+para.Wcortex);

			// without Golgi cell feedback
			if(para.gocmodeflag==0){
				if(para.actfuncflag==1 && para.actfuncflag_GC==2){
					// linear RNN & cereb
					para.Wrank1=para.Wrank1*arma::diagmat(para.gain_vector)*para.Wctxgc;
					para.Weff=para.WrankN+para.Wrank1;
				}
				else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					para.Wctxgc=para.Wctxgc;
				}
			}
			// with Golgi cell feedback
			else if (para.gocmodeflag==1){
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
				}
			}
			// with feedforward Golgi cell inhibition
			else if (para.gocmodeflag==3)	{
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					para.WrankN=(-para.I+para.Wcortex);
				}
			}

			nnet = static_cast<size_t>(para.Ncx);
			nn2 = static_cast<size_t>(para.Ncx+para.Ncx);
			FLOQsys = {RK_Floquet_cortex_cereb_rhs_vector, NULL, nn2, &para};
		}
		// *************************************
		// non-instantaneous GCs
		// *************************************
		else{
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
			}

			nnet = static_cast<size_t>(para.Ncx+para.Ngc);
			nn2 = static_cast<size_t>(para.Ncx+para.Ngc+para.Ncx+para.Ngc);
			FLOQsys = {RK_Floquet_cortex_cereb_dyn_rhs_vector, NULL, nn2, &para};
		}
	}

	// *****************************
	// integrate system
	// *****************************
	gsl_odeiv2_driver * FLOQdriver = gsl_odeiv2_driver_alloc_y_new (&FLOQsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs/10.0,para.INT_delta_rel/10.0);

	std::vector<double> Xnet(nnet);
	std::vector<double> X(nn2);
	int status;

	// random perturbation vector
	vm = arma::zeros<arma::vec>(nnet);
	for(int i= 0; i < nnet; i++) {
		vm(i)=gsl_ran_gaussian_ziggurat (para.r, 1.0);
	}

	// Normalise initial perturbation to unit norm
	vm = vm/arma::norm(vm);

	// assign state variables
	if(para.taugcmodeflag==0){
		for(int i= 0; i < para.Ncx; i++) {
			X[i]=XIN[i];
			X[para.Ncx+i]=vm(i);
		}
	}
	else if(para.taugcmodeflag>0){
		for(int i= 0; i < para.Ncx+para.Ngc; i++) {
			X[i]=XIN[i];
			X[para.Ncx+para.Ngc+i]=vm(i);
		}
	}

	double log_growth = 0.0;
	double time = 0.0;

	para.closedloop=para.floquetmodeflag;	// closed/open loop according to chosen parameter
	para.testphaseflag=0;					// no noise
	para.RKflag=1;							// Runke-Kutta integration
	log_file << "\t \t \t FM computation phase" << "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	for (int k = 0; k < para.N_FMperiods; k++) {
		// integrate one period
		gsl_odeiv2_driver_apply(FLOQdriver, &time, time + para.Tmax_FM*para.dt, X.data());

		// final state vector
		for(int i= 0; i < nnet; i++) {
			Xnet[i]=X[i];							// final state vector
		}

		// final perturbation vector
		arma::vec vT(X.data() + nnet, nnet);

		// compute tangent vector at current time
		std::vector<double> dX_dt(nnet);
		if(para.algswitch==1) {
			RK_cortex_rhs(time, Xnet.data(), dX_dt.data(), &para);
		}
		else {
			if(para.taugcmodeflag==0){
				RK_cortex_cereb_rhs(time, Xnet.data(), dX_dt.data(), &para);
			}
			else if(para.taugcmodeflag>0){
				RK_cortex_cereb_dyn_rhs(time, Xnet.data(), dX_dt.data(), &para);
			}
		}
		arma::vec tangent_vec(dX_dt.data(), nnet);

		double tangent_norm2 = arma::dot(tangent_vec, tangent_vec);
		// project out trivial direction
		arma::vec v_perp = vT - arma::dot(vT, tangent_vec) / arma::dot(tangent_vec, tangent_vec) * tangent_vec;

		// accumulate log growth
		log_growth += std::log(arma::norm(v_perp));

		// renormalise for next period
		v_perp /= arma::norm(v_perp);
		std::copy(v_perp.begin(), v_perp.end(), X.data() + nnet);
	}
	log_file << "\t \t \t End of FM computation phase, T =  "<< para.N_FMperiods*(para.Tmax_FM)*para.dt << "\t Ttot =  "<< time << std::endl;

	double mu_max = std::exp(log_growth / para.N_FMperiods);
	double LE_max = log_growth / (para.N_FMperiods*para.Tmax_FM*para.dt);

	if(para.learnflag<99){
		std::cout << "Largest non-trivial Floquet multiplier:" << std::endl;
		std::cout << mu_max << std::endl;
		std::cout << "Largest finite time lyapunov exponent:" << std::endl;
		std::cout << LE_max << std::endl;
		std::cout << std::endl;
	}

	if(para.algswitch==1)	json_measures["Max_Floquet_multiplier_ctx"].push_back(mu_max);
	else					json_measures["Max_Floquet_multiplier_cereb"].push_back(mu_max);

	gsl_odeiv2_driver_free (FLOQdriver);

}
else{
	if(para.algswitch==1)	json_measures["Max_Floquet_multiplier_ctx"].push_back(NAN);
	else					json_measures["Max_Floquet_multiplier_cereb"].push_back(NAN);
}

	return;
}
