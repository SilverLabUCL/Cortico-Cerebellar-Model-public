void set_taus(netpara &para, double dt) {
	//*******************************
	// RNN taus
	//*******************************
	para.sqt=sqrt(dt);
	para.dttau_distrb=dt/para.tau_distrb;
	para.IIdt_distb=(1.0-para.dttau_distrb);
	para.sqDT=para.sqt/sqrt(para.tau_distrb);

	//*******************************
	// GC taus
	//*******************************
	para.dttauGC_distb=dt/para.tauGC_distrb;
	para.IIdtGC_distb=(1.0-para.dttauGC_distb);
	//~ para.sqtauGC_distb=sqrt(para.tauGC_distrb);
	para.sqDTGC=para.sqt/sqrt(para.tauGC_distrb);

	return;
}
