void update_rates_ctx(double t, const double X[], netpara &para) {

	// take system state from integrator and calculate firing rates
	for(int i= 0; i < para.Ncx; i++) {
		para.nuCx[i]=(*para.Ptr_ctxIO)(X[i],1.0);
	}

	return;
}
