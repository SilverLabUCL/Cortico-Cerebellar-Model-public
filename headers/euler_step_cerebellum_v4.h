void euler_step_cerebellum(arma::colvec& nuCx, arma::colvec& GC, arma::colvec& hGC, double &time, netpara &para) {

	double GoC=0;
	double sigNsq2=para.sigmaNgc/sqrt(2.0);
	arma::vec MF(para.Nmf);
	arma::vec input_intf(para.Ncomp_intf);

	// MFs = sample from cortical activity
	for (int j=0; j<para.Nmf; j++) 	MF[j]=nuCx[para.Ncx_to_MF_idx[j]];

	if(para.gocmodeflag==1 || para.gocmodeflag==2){
		for (int i=0; i<para.Ngc; i++) {
			GoC+=GC[i];
		}
		GoC=GoC/double(para.Ngc);
	}
	else if(para.gocmodeflag==3){
		GoC=arma::stddev(MF);
	}
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			//~ input_intf[j]=func_target_intf(t, j, para);
			input_intf[j]=func_target_intf(time, j, para);
		}
	}

	// evolve GCs
	for (int i=0; i<para.Ngc; i++) {
		//	pre-synaptic inputs of GC i
		double inputs=0;
		for (int j=0; j<para.Nsyn; j++) {
			//~ inputs+=MF[Wmf[i].idx[j]];
			inputs+=MF[para.Wmf[i].idx[j]];
			//~ inputs+=Wmf[i].weight[j]*MF[Wmf[i].idx[j]];
		}
		inputs=inputs/para.sqNsyn;
		if(para.gocmodeflag>0) inputs+=-para.Wgoc_distrb[i]*GoC;

		if(para.taugcmodeflag==0){
			// instananeous GCs
			hGC[i]=inputs;
			// constant noise
			//~ if(para.noiseflag==0 && para.testphaseflag==1 && para.sigmaNgc>0)	hGC[i]=hGC[i] + sigNsq2*gsl_ran_gaussian_ziggurat(para.r,1.0);
		}
		else if(para.taugcmodeflag>0){
			// GCs with identical or distributed tauGC
			// hGC[i]= (1-dt/tauGC)*hGC[i]+dt/tauGC*inputs;    //+ para.sqDT*sigmaN*xi;
			hGC[i]= para.IIdtGC_distb[i]*hGC[i]+para.dttauGC_distb[i]*inputs;
			// constant noise
			//~ if(para.noiseflag==0 && para.testphaseflag==1 && para.sigmaNgc>0)	hGC[i]=hGC[i] + para.sqDTGC[i]*gsl_ran_gaussian_ziggurat(para.r,para.sigmaNgc);
		}
	}

	// add noise
	if(para.noiseflag==0 && para.testphaseflag==1 && para.sigmaNgc>0){
		// instantaneous GCs
		if(para.taugcmodeflag==0){
			for (int i=0; i<para.Ngc; i++) {
				// constant noise
				hGC[i]=hGC[i] + sigNsq2*gsl_ran_gaussian_ziggurat(para.r,1.0);
			}
		}
		// ODE GCs with identical or distributed tauGC
		else if(para.taugcmodeflag>0){
			for (int i=0; i<para.Ngc; i++) {
				// constant noise
				hGC[i]=hGC[i] + para.sqDTGC[i]*gsl_ran_gaussian_ziggurat(para.r,para.sigmaNgc);
			}
		}

	}

	// GC rates
	for (int i=0; i<para.Ngc; i++) {
		GC[i]=para.gain_vector[i]*(*para.Ptr_gcIO) (hGC[i]-para.threshold_vector[i],para.sp_alpha);

		// add interference signal
		if(para.intrfflag==2){
			//~ for (int j=0; j<para.Ncomp_intf; j++)	GC[i]=GC[i]+interfsignal(j,t)*para.Intf_weights(i,j);
			for (int j=0; j<para.Ncomp_intf; j++)	GC[i]=GC[i]+input_intf[j]*para.Intf_weights.at(i,j);
		}

	}

	// GC filtering
	if(para.filterflag>0)	GC=filter_GCs(GC, para);

	return;
}
