#ifndef ALLSTRUCTS_H
#define ALLSTRUCTS_H

struct presyn {
	int Nsyn;
	arma::uvec idx;
	arma::vec weight;

    // Full constructor
	presyn(const arma::uvec& idx_, const arma::vec& weight_)
		: Nsyn(idx_.n_elem), idx(idx_), weight(weight_) {}

	// Weight-only constructor
	presyn(const arma::vec& weight_)
		: Nsyn(weight_.n_elem), idx(), weight(weight_) {}
};

struct RLS {
	int N;
	arma::vec dwout;
	arma::mat P;
	double LR;
};

struct momenta {
	double f;
	double hctx_mu,hctx_sig;
	double hgc_mu,hgc_sig;
	double Ctx_mu,Ctx_sig;
	double GC_mu, GC_sig;
};

struct netpara {
	int Ncx,Ngc,Nmf,Nsyn,Ncomp,Ncomp_noise,Ncomp_noise_RNN,Ncomp_intf,Nwrite,Npara,N_FMperiods,Nreal,Nrep,seed;
	int Tmax_learn,Tmax_test,Tmax_test_euler,Tmax_FT,Tmax_trans,Tmax_FM,Tmax_Fadj;
	int algswitch,RLSlearn_every,write_every,closedloop;

	double sigmaN,sigmaNgc,sqNsyn,alpha,Wgoc,Wgoc_SD,pcon,g,sigmaLS,FACsigmaLS_ctx,FACsigmaLS_gc,tol_GCgain,GCgainmax,GCmax;
	double sqt,delta_eta,tau,delta_r2_reg,tau_avrg,tau_std,tauRNNcutoff,tauGC_avrg,tauGC_std,tauGCcutoff,dt,dt_euler,DT,DT_euler,ftarget,GCsig_target,GCavrg_target,GCgain,GCthreshold,GCgain_base,GCthreshold_base,fc_avrg,fc_SD,rho_Wfb12,rho_Wfb13;
	double sigAmp,sigAmp_intf,lambda,phi1,ampmod,sp_alpha,FCcomp_tol;
	double Z,PC;

	double duration_learn,duration_test,duration_trans,duration_FT,signalavrg,signalvar,error_signed_train;
	double INT_delta_abs,INT_delta_rel;

	// flags
	int actfuncflag,actfuncflag_GC,fcodingflag,floquetflag,floquetmodeflag,gainflag,gocmodeflag,intrfflag,learnflag,loop_learning_flag,loop_test_flag,LRflag,noiseflag,RKflag,ROwriteflag,signalflag,simplemapflag,taudistrbflag,taugcmodeflag,testphaseflag;
	int cortexonlyflag;

	// strings
	std::string activation_function,activation_function_GC,fcodingmode,floquet,floquetmode,gainmode,GoC_mode,interference_signal,learningrule,learnmode,loop_test,MF_GC_connect,noisescaling,path,path_connect,ROwrite,signaltype,tau_mode,taugc_mode;
	std::string folder,version;

	//arrays
	arma::vec nuCx;
	arma::vec nuGC;
	arma::vec hGC;

	arma::vec tau_distrb;
	arma::vec inv_tau_distrb;
	arma::vec IIdt_distb;
	arma::vec dttau_distrb;
	arma::vec tauGC_distrb;
	arma::vec inv_tauGC_distrb;
	arma::vec IIdtGC_distb;
	arma::vec dttauGC_distb;
	arma::vec sqDT;
	arma::vec sqDTGC;
	arma::vec inv_tauALL_distrb;

	arma::vec freqs;
	arma::vec amps;
	arma::vec phases;
	arma::vec freqs_intf;
	arma::vec amps_intf;
	arma::vec phases_intf;

	arma::vec freqs_func;
	arma::vec freqs_intf_func;
	arma::vec amps_func;
	arma::vec phases_func;
	arma::vec amps_intf_func;
	arma::vec phases_intf_func;

	arma::vec fc_target;
	arma::vec gain_vector;
	arma::vec threshold_vector;

	arma::uvec Ncx_to_MF_idx;		// cx-to-MF index

	arma::vec Wgoc_distrb;
	arma::vec Wgoc_distrb_scaled;

	arma::mat rho_Wfb;
	arma::mat Intf_weights;

	std::vector<presyn> Wcx;		// cx-to-cx
	std::vector<presyn> Wmf;		// cx-to-gc
	arma::vec Wfb;					// output-unit-to-cx
	arma::mat Wfb_interf;			// interference-to-cx
	arma::vec Wout;					// cx-to-output-unit
	arma::vec Wout_PC;				// cx-to-output-unit

	arma::mat Weff;
	arma::mat Wcortex;
	arma::mat WrankN;
	arma::mat Wrank1;
	arma::mat Wctxgc;
	arma::mat I;
	arma::mat Igc;
	arma::mat CD;

	//	function pointer
	double (*Ptr_ctxIO)(double,double);
	double (*Ptr_gcIO)(double,double);

	// randum numbers
	gsl_rng * r;

};

struct FTcomp {
	arma::vec a0;
	arma::mat an;
	arma::mat bn;
	arma::mat freq_comp;
	arma::mat S_signal;
	arma::mat S_noise;
	arma::mat S_noise_noa0;
	arma::mat S_total;
	arma::mat S_signal_centered;
	arma::mat S_noise_centered;
	arma::mat S_noise_noa0_centered;
	arma::mat S_total_centered;
	arma::mat S_signal_sqN;
	arma::mat S_noise_sqN;
	arma::mat S_total_sqN;

	arma::mat Angles_signal;
	arma::mat Angles_all;

	arma::vec nn_signal;
	arma::vec nn_noise;
	arma::vec nn_total;

	arma::vec subangles;

	double dim_signal;
	double cond_signal,cond_total;
	double sineprod_signal;
	int Ms,Mn;
};

#endif
