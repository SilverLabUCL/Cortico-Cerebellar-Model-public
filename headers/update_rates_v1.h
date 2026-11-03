void update_rates(double t, const double X[], netpara &para) {

	// take system state from integrator and calculate firing rates

	double GoC;
	arma::vec input_intf(para.Ncomp_intf);
	arma::vec MF(para.Nmf);

	// RNN firing rates
	for(int i= 0; i < para.Ncx; i++) {
		para.nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	// GC inputs
	if(para.taugcmodeflag==0){
		// MFs = sample from cortical activity
		for (int j=0; j<para.Nmf; j++) 	MF[j]=para.nuCx[para.Ncx_to_MF_idx[j]];

		// Golgi cell
		if(para.gocmodeflag==1 || para.gocmodeflag==2){
			GoC=arma::mean(para.nuGC);
		}
		else if(para.gocmodeflag==3){
			GoC=arma::stddev(MF);
		}

		// evolve GCs
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
			para.hGC[i]=inputs;
		}

	}
	else if(para.taugcmodeflag>0){
		for (int i=0; i<para.Ngc; i++) {
			para.hGC[i]=X[i+para.Ncx];
		}
	}

	// GC interference signal
	if(para.intrfflag==2){
		for (int j=0; j<para.Ncomp_intf; j++){
			input_intf[j]=func_target_intf(t, j, para);
		}
	}

	// GC rates
	for (int i=0; i<para.Ngc; i++) {
		para.nuGC[i]=para.gain_vector[i]*(*para.Ptr_gcIO) (para.hGC[i]-para.threshold_vector[i],para.sp_alpha);

		// add GC interference signal
		if(para.intrfflag==2){
			for (int j=0; j<para.Ncomp_intf; j++)	para.nuGC[i]=para.nuGC[i]+input_intf[j]*para.Intf_weights.at(i,j);
		}
	}

	// GC filtering
	if(para.filterflag>0)	para.nuGC=filter_GCs(para.nuGC, para);


	return;
}
