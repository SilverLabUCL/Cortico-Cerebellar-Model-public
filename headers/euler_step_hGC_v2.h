void euler_step_hGC(arma::colvec& nuCx, double hGC[], arma::colvec& MF, netpara &para) {

	double a;
	double GoC;
	// MFs = sample from cortical activity
	for (int j=0; j<para.Nmf; j++) {
		MF[j]=nuCx[para.Ncx_to_MF_idx[j]];
	}
	if(para.gocmodeflag==3)		GoC=arma::stddev(MF);

	// evolve GCs
	for (int i=0; i<para.Ngc; i++) {
		//	pre-synaptic inputs of GC i
		a=0;
		for (int j=0; j<para.Nsyn; j++) {
			a+=MF[para.Wmf[i].idx[j]];
			//~ a+=Wmf[i].weight[j]*MF[Wmf[i].idx[j]];
		}

		if(para.taugcmodeflag==0){
			// instananeous GCs
			if(para.gocmodeflag==3)		hGC[i]=a/para.sqNsyn - para.Wgoc*GoC;
			else 						hGC[i]=a/para.sqNsyn;
		}
		else if(para.taugcmodeflag>0){
			// GCs with identical or distributed tauGC
			// hGC[i]= (1-dt/tauGC)*hGC[i]+dt/tauGC*a/para.sqNsyn;
			if(para.gocmodeflag==3)		hGC[i]=para.IIdtGC_distb[i]*hGC[i]+para.dttauGC_distb[i]*(a/para.sqNsyn - para.Wgoc*GoC);
			else						hGC[i]=para.IIdtGC_distb[i]*hGC[i]+para.dttauGC_distb[i]*a/para.sqNsyn;
		}
	}

	return;
}
