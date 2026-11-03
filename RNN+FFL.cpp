// added boost::program_options

#include <vector>
#include <armadillo>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_odeiv2.h>
#include <sys/stat.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "headers/json.hpp"
#include "headers/functions.h"
#include "headers/allstructs_v5.h"

#include "headers/target_signal_v6.h"
#include "headers/euler_updates_v1.h"
#include "headers/RK_rhs_v1.h"
#include "headers/RK_weights_rhs_v1.h"
#include "headers/update_rates_v1.h"
#include "headers/RLS_step_cerebellum_allGCs_v3.h"
#include "headers/calc_FT_v7.h"
#include "headers/calc_floquet_RK_vector_v1.h"
#include "headers/calc_artificial_overtones_v4.h"
#include "headers/LS_solution_v4.h"
#include "headers/adjust_thresholds_gains_v2.h"
#include "headers/generate_weights_v3.h"
#include "headers/output_calculations_v2.h"
#include "headers/calc_linCxVar.h"
#include "headers/set_scan_variables.h"
#include "headers/set_parameters_flags_v2.h"
#include "headers/write_parameters_v2.h"
#include "headers/write_measures_final_v4.h"
#include "headers/set_taus_v1.h"
#include "headers/write_single_v2.h"

int main (int argc, char* argv[]) {

//------------------------------------------------------------------------
//	set parameters, flags, write to file & print some output
//------------------------------------------------------------------------
	struct netpara para;
	para.version="RNN+FFL_public_v1.0";

	// dummy parameters
	para.cortexonlyflag=0;

	int rtidx=set_parameters_flags(para,argc,argv);
	if(rtidx==1) return(0);

//------------------------------------------------------------------------
// 	set up random number generator
//------------------------------------------------------------------------
	srand ( time(NULL) );

  	const gsl_rng_type * T;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	para.r = gsl_rng_alloc (T);

	gsl_rng_env_setup();
  	gsl_rng_set(para.r,para.seed); 			// seed the random number generator
//------------------------------------------------------------------------
//	set up parameter array for scan
//------------------------------------------------------------------------
	arma::vec SCANPARA = set_scan_variables(para.Npara, para.learnflag, para);
//------------------------------------------------------------------------
//	create json object for parameter output & single trial output
//------------------------------------------------------------------------
	nlohmann::json json_parameters;
	nlohmann::json json_trial;
//------------------------------------------------------------------------
// 	prepare output files
//------------------------------------------------------------------------
	for(int np= 0; np < para.Npara; np++) {
		para.folder=para.path+"parameter_"+ std::to_string (np+1);
		struct stat sb;
		if (stat(para.folder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
			if(para.Npara<20) std::cout << "[Writing in existing folder]" << "\t" << para.folder << std::endl;
		}
		else{
			std::string dircommand=std::string("mkdir ") + para.folder;
			int systemRet = system(dircommand.c_str());
			if(systemRet == -1){
				std::cerr <<"Error: Folder could not be created."<< std::endl;
				exit(EXIT_FAILURE);
			}
			else if(para.Npara<20) std::cout << "[Creating folder]"  << "\t" << para.folder << std::endl;
		}
	}
	if(para.Npara>=20){
		std::cout << "[Creating/Writing into " << para.Npara << " folders]" << std::endl;
	}
	std::cout << std::endl;

	std::ofstream nuCx_learn_file(para.path+"nuCx_learn.dat");		if (!nuCx_learn_file) exit(EXIT_FAILURE);
	std::ofstream nuCx_test_file(para.path+"nuCx_test.dat");		if (!nuCx_test_file) exit(EXIT_FAILURE);
	std::ofstream nuCx_FT_file(para.path+"nuCx_FT.dat");			if (!nuCx_FT_file) exit(EXIT_FAILURE);

	std::ofstream nuGC_learn_file(para.path+"GC_learn.dat");			if (!nuGC_learn_file) exit(EXIT_FAILURE);
	std::ofstream nuGC_test_file(para.path+"GC_test.dat");			if (!nuGC_test_file) exit(EXIT_FAILURE);
	std::ofstream nuGC_FT_file(para.path+"GC_FT.dat");				if (!nuGC_FT_file) exit(EXIT_FAILURE);
	std::ofstream hGC_FT_file(para.path+"hGC_FT.dat");				if (!hGC_FT_file) exit(EXIT_FAILURE);
	std::ofstream Inputs_file(para.path+"Inputs.dat");				if (!Inputs_file) exit(EXIT_FAILURE);

	std::ofstream freq_file(para.path+"frequency.dat");			if (!freq_file) exit(EXIT_FAILURE);
	std::ofstream amp_file(para.path+"amp.dat");				if (!amp_file) exit(EXIT_FAILURE);
	std::ofstream phase_file(para.path+"phase.dat");			if (!phase_file) exit(EXIT_FAILURE);
	std::ofstream freq_intf_file(para.path+"frequency_intf.dat");		if (!freq_intf_file) exit(EXIT_FAILURE);
	std::ofstream amp_intf_file(para.path+"amp_intf.dat");			if (!amp_intf_file) exit(EXIT_FAILURE);
	std::ofstream phase_intf_file(para.path+"phase_intf.dat");		if (!phase_intf_file) exit(EXIT_FAILURE);

	std::ofstream angles_ctx_file(para.path+"angles_ctx.dat");		if (!angles_ctx_file) exit(EXIT_FAILURE);
	std::ofstream angles_cereb_file(para.path+"angles_cereb.dat");	if (!angles_cereb_file) exit(EXIT_FAILURE);
	std::ofstream S_total_ctx_file(para.path+"stotal_ctx.dat");		if (!S_total_ctx_file) exit(EXIT_FAILURE);
	std::ofstream S_total_cereb_file(para.path+"stotal_cereb.dat");	if (!S_total_cereb_file) exit(EXIT_FAILURE);

	std::ofstream error_signed_train_file(para.path+"error_signed_train.dat");	if (!error_signed_train_file) exit(EXIT_FAILURE);
	std::ofstream error_full_train_file(para.path+"error_full_train.dat");		if (!error_full_train_file) exit(EXIT_FAILURE);
	std::ofstream error_signed_test_file(para.path+"error_test.dat");			if (!error_signed_test_file) exit(EXIT_FAILURE);
	std::ofstream LR_file(para.path+"LR.dat");									if (!LR_file) exit(EXIT_FAILURE);

	std::ofstream simcount_file(para.path+"simcount.dat");				if (!simcount_file) exit(EXIT_FAILURE);

	std::ofstream single_ctxunit_files(para.path+"single_ctxunit.dat");	if (!single_ctxunit_files) exit(EXIT_FAILURE);
	std::ofstream single_gc_files(para.path+"single_gc.dat");			if (!single_gc_files) exit(EXIT_FAILURE);
	std::ofstream single_output_files(para.path+"single_output.dat");	if (!single_output_files) exit(EXIT_FAILURE);
	std::ofstream single_ctxunit_FT_files(para.path+"single_ctxunit_FT.dat");	if (!single_ctxunit_FT_files) exit(EXIT_FAILURE);
	std::ofstream single_gc_FT_files(para.path+"single_gc_FT.dat");			if (!single_gc_FT_files) exit(EXIT_FAILURE);
	std::ofstream single_output_FT_files(para.path+"single_output_FT.dat");	if (!single_output_FT_files) exit(EXIT_FAILURE);

	std::ofstream alg_steps_file(para.path+"alg_steps.dat");			if (!alg_steps_file) exit(EXIT_FAILURE);

	std::ofstream PC_test_files[para.Npara];

	for(int np= 0; np < para.Npara; np++) {
		para.folder=para.path+"parameter_"+ std::to_string (np+1) + "/";

		PC_test_files[np].open(para.folder+"PC_test.dat", std::ofstream::out | std::ofstream::trunc);
		PC_test_files[np].close();
	}

//	**************************************************
//	Npara loop
//	**************************************************
	for(int np= 0; np < para.Npara; np++) {
		alg_steps_file << "Npara = " << np+1 << std::endl;

		//------------------------------------------------------------------------
		//	create json objects for output of measures, stats & activity traces
		//------------------------------------------------------------------------
		nlohmann::json json_measures;
		nlohmann::json json_stats;

		// open output files
		para.folder=para.path+"parameter_"+ std::to_string (np+1) + "/";

		PC_test_files[np].open(para.folder+"PC_test.dat",std::ios::app);

		// set parameters
		if(para.learnflag==100){
			para.sigmaN=SCANPARA(np);
			para.sigmaNgc=para.sigmaN;
			std::cout << "sigmaN = " << para.sigmaN  << "\t sigmaNgc = " << para.sigmaNgc<< std::endl;
		}
		else if(para.learnflag==101){
			para.GCthreshold_base=SCANPARA(np);
			std::cout << "GC threshold = " << para.GCthreshold_base << std::endl;
		}
		else if(para.learnflag==102){
			para.Ncx=SCANPARA(np);
			para.Nmf=para.Ncx;
			std::cout << "Ncx = " << para.Ncx << std::endl;
		}
		else if(para.learnflag==103){
			para.Ngc=SCANPARA(np);
			std::cout << "Ngc = " << para.Ngc << std::endl;
		}
		else if(para.learnflag==105){
			para.sigmaNgc=SCANPARA(np);
			std::cout << "sigmaNgc = " << para.sigmaNgc << std::endl;
		}
		else if(para.learnflag==106){
			para.Wgoc=SCANPARA(np);
			std::cout << "Wgoc = " << para.Wgoc << std::endl;
		}
		else if(para.learnflag==107){
			para.Ncomp=SCANPARA(np);
			std::cout << "Ncomp = " << para.Ncomp << std::endl;
		}
		else if(para.learnflag==108){
			para.sigAmp=SCANPARA(np);
			std::cout << "sigAmp = " << para.sigAmp << std::endl;
		}
		else if(para.learnflag==109){
			para.lambda=SCANPARA(np);
			std::cout << "lambda = " << para.lambda << std::endl;
		}
		else if(para.learnflag==110){
			para.duration_learn=SCANPARA(np);
			std::cout << "Tlearn = " << para.duration_learn << std::endl;
		}
		else if(para.learnflag==113){
			para.g=SCANPARA(np);
			std::cout << "g = " << para.g << std::endl;
		}
		else if(para.learnflag==114){
			para.delta_eta=SCANPARA(np);
			std::cout << "delta_eta = " << para.delta_eta << std::endl;
		}
		else if(para.learnflag==118){
			para.phi1=SCANPARA(np);
			std::cout << "Phi1 = " << para.phi1 << std::endl;
		}
		else if(para.learnflag==119){
			para.ampmod=SCANPARA(np);
			std::cout << "Ampmod = " << para.ampmod << std::endl;
		}
		else if(para.learnflag==120){
			para.sigmaLS=SCANPARA(np);
			std::cout << "sigmaLS = " << para.sigmaLS << std::endl;
		}
		else if(para.learnflag==123){
			para.fc_avrg=SCANPARA(np);
			std::cout << "f_coding average = " << para.fc_avrg << ", f_coding SD = "<< para.fc_SD << std::endl;
		}
		else if(para.learnflag==124){
			para.sp_alpha=SCANPARA(np);
			std::cout << "Softplus sharpness alpha = " << para.sp_alpha << std::endl;
		}
		else if(para.learnflag==128){
			para.rho_Wfb12=SCANPARA(np);
			std::cout << "rho_Wfb_12 = " << para.rho_Wfb12 << std::endl;
		}
		else if(para.learnflag==129){
			para.sigAmp_intf=SCANPARA(np);
			std::cout << "sigAmp_intf = " << para.sigAmp_intf << std::endl;
		}
		else if(para.learnflag==131){
			para.GCgain_base=SCANPARA(np);
			std::cout << "GCgain_base = " << para.GCgain_base << std::endl;
		}
		else if(para.learnflag==132){
			para.rho_Wfb13=SCANPARA(np);
			std::cout << "rho_Wfb_13 = " << para.rho_Wfb13 << std::endl;
		}
		else if(para.learnflag==134){
			para.GCgainmax=SCANPARA(np);
			std::cout << "GCgainmax = " << para.GCgainmax << std::endl;
		}

	para.dt=para.DT/para.lambda;
	para.dt_euler=para.DT_euler/para.lambda;
	para.Tmax_learn = round(para.duration_learn/para.DT);
	para.Tmax_test = round(para.duration_test/para.DT);
	para.Tmax_test_euler = round(para.duration_test/para.DT_euler);
	para.Tmax_FT = round(para.duration_FT/para.DT);
	para.Tmax_trans = round(para.duration_trans/para.DT);
	para.Tmax_FM = round(1.0/para.DT);						// one period for FM input arrays
	para.Tmax_Fadj = round(2.0/para.DT);					// two periods

	int status;

	double hcx[para.Ncx],hcx_save[para.Ncx],hcxini[para.Ncx];
	double X[para.Ncx+para.Ngc],X_save[para.Ncx+para.Ngc];
	double time,error_signed_test,RMSE_train,RMSE_test;

	struct RLS RLS_cereb;
	RLS_cereb.N=para.Ngc;
	RLS_cereb.dwout=arma::zeros<arma::vec>(para.Ngc);
	RLS_cereb.P=arma::zeros<arma::mat>(para.Ngc,para.Ngc);

	arma::mat RATES_FT(para.Ncx, para.Tmax_FT);
	arma::mat RATES_FTtrans(para.Ncx, para.Tmax_trans);
	arma::mat RATES_FT_GC(para.Ngc, para.Tmax_FT);
	arma::mat RATES_LS_GC(para.Ngc, para.Tmax_learn);
	arma::mat INPUTS_Fadj_GC(para.Ngc, para.Tmax_Fadj);
	arma::mat RATES_Fadj(para.Ncx, para.Tmax_Fadj);

	int Ninput;
	if(para.taugcmodeflag==0){
		Ninput=para.Ncx;
	}
	else if(para.taugcmodeflag>0){
		Ninput=para.Ncx+para.Ngc;
	}
	arma::mat INPUTS_OPEN(Ninput, para.Tmax_FM);
	arma::mat INPUTS_CLOSED(Ninput, para.Tmax_FM);

	arma::colvec fsparse(para.Ngc);
	arma::colvec fsparse_pop(para.Tmax_FT);
	arma::colvec nuGC_avrg_pop(para.Tmax_FT);
	arma::colvec OUTvec(para.Tmax_learn);
	arma::colvec Target(para.Tmax_learn);

	arma::running_stat<double> stats_hcx;
	arma::running_stat<double> stats_nuCx;
	arma::running_stat<double> stats_hgc;
	arma::running_stat<double> stats_GC;
	arma::running_stat<double> stats_rmse;
//------------------------------------------------------------------------
// 	generate target signal(s)
//------------------------------------------------------------------------
	alg_steps_file << "\t \t write target signal(s) to: "<< para.path+"parameter_"+ std::to_string (np+1) + "/" << std::endl;
	target_signal(para,freq_file,amp_file,phase_file,freq_intf_file,amp_intf_file,phase_intf_file,alg_steps_file);

//------------------------------------------------------------------------
// 	find population activity momenta of naive network
//------------------------------------------------------------------------
	alg_steps_file << "\t find population activity momenta of naive network" << std::endl;
	para.GCthreshold=para.GCthreshold_base/para.sqNsyn;

	// if RNN is linear, adjust thresholds in function of population activity variance
	double CxVar;
	if(para.actfuncflag==1){
		CxVar=calc_linCxVar(para);
		para.GCthreshold=para.GCthreshold*sqrt(CxVar);
	}

//	*************************
// 	Set sigmaLS multiplier
//	*************************
	para.FACsigmaLS_ctx=para.sigAmp;
	para.FACsigmaLS_gc=para.sigAmp;
	if(para.gainflag==4){
		para.FACsigmaLS_gc=para.sigAmp*para.GCgain;
	}

//------------------------------------------------------------------------
// 	write parameters to file
//------------------------------------------------------------------------
	write_parameters(json_parameters,para);

//	**************************************************
//	Nreal loop
//	**************************************************
	for (int rr = 0; rr < para.Nreal; rr++){
		if( ( (rr+1) % 5) ==0) std::cout << "Realisation number " << rr+1 << std::endl;
		alg_steps_file << "\t Nreal = " << rr+1 << std::endl;

//	--------------------------------------------------
//	set up connectivity
//	--------------------------------------------------
	generate_weights(para);
	alg_steps_file << "\t \t generate new weights"<< std::endl;
//	--------------------------------------------------
// initial cortical state
//	--------------------------------------------------
	for(int i= 0; i < para.Ncx; i++) {
		hcxini[i]=0.5*gsl_ran_gaussian_ziggurat(para.r,1);
		//~ hcxini[i]=0;
	}
	if (para.learnflag<99) {
		for(int i= 0; i < para.Ncx; i++) json_trial["hcx_inistate"].push_back(hcxini[i]);
	}
	alg_steps_file << "\t \t generate initial conditions"<< std::endl;

//------------------------------------------------------------------------
// 	begin integration of the system
//------------------------------------------------------------------------

	//	*************************
	// 	set up GSL integrator
	//	*************************
	size_t nrnn = static_cast<size_t>(para.Ncx);
	size_t nrnngc = static_cast<size_t>(para.Ncx+para.Ngc);
	// RNN only
  	gsl_odeiv2_system RNNsys = {RK_cortex_rhs, NULL, nrnn, &para};

  	// RNN + cereb
  	gsl_odeiv2_system RNNGCsys;
	if(para.taugcmodeflag==0){
		RNNGCsys= {RK_cortex_cereb_rhs, NULL, nrnn, &para};
	}
	else if(para.taugcmodeflag>0){
		RNNGCsys= {RK_cortex_cereb_dyn_rhs, NULL, nrnngc, &para};
	}

	gsl_odeiv2_driver * RNNdriver = gsl_odeiv2_driver_alloc_y_new (&RNNsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

	gsl_odeiv2_driver * RNNGCdriver = gsl_odeiv2_driver_alloc_y_new (&RNNGCsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

	//	*************************
	// 	set initial conditions
	//	*************************
	para.nuCx= arma::zeros<arma::vec>(para.Ncx);
	para.hGC= arma::zeros<arma::vec>(para.Ngc);
	para.nuGC= arma::zeros<arma::vec>(para.Ngc);
	para.Wout_PC= arma::zeros<arma::vec>(para.Ngc);
	for(int i= 0; i < para.Ncx; i++) {
		hcx[i]=hcxini[i];
	}
	// set GC thresholds and Gains
	para.gain_vector=arma::ones<arma::vec>(para.Ngc)*para.GCgain;
	para.threshold_vector=arma::ones<arma::vec>(para.Ngc)*para.GCthreshold;
	para.Intf_weights=arma::zeros<arma::mat>(para.Ncx,para.Ncomp_intf);
	time=0;
	int tt=0;

	//	*************************
	// 	initial transient phase (always open loop / driven)
	//	*************************
	// set loop parameter
	para.closedloop=0;		// open loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t first transient phase"<< "\t RKflag = "<< para.RKflag << " \t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;

	tt=para.Tmax_trans;
	status = gsl_odeiv2_driver_apply (RNNdriver, &time, (para.Tmax_trans-1)*para.dt, hcx);
	if (status != GSL_SUCCESS){
		printf ("Error in first transient phase, return value=%d\n", status);
		break;
	}
	alg_steps_file << "\t \t \t End of first transient phase, T =  "<< para.Tmax_trans*para.dt  << "\t Ttot =  "<< time << std::endl;

	for(int i= 0; i < para.Ncx; i++) {
		X[i]=hcx[i];
	}
	for(int i= 0; i < para.Ngc; i++) {
		X[i+para.Ncx]=0.0;
	}

	//	*************************
	// f coding adjustment phase (always open loop / driven)
	//	*************************
	// set loop parameter
	para.closedloop=0;		// open loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t Coding adjustment phase"<< "\t RKflag = "<< para.RKflag << "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	for(int t= 0; t < para.Tmax_Fadj; t++){															// Tmax_Fadj showld not be larger than Tmax_trans!!!!

		if(para.taugcmodeflag==0){
			status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, tt*para.dt, hcx);
			update_rates(time, hcx, para);
		}
		else if(para.taugcmodeflag>0){
			status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, tt*para.dt, X);
			for(int i= 0; i < para.Ngc; i++) {
				para.hGC[i]=X[i+para.Ncx];
			}
			update_rates(time, X, para);
		}

		if (status != GSL_SUCCESS){
			printf ("Error in coding adjustment phase, return value=%d\n", status);
			break;
		}

		// set activity of output equal to target to ensure continuity with subsequent phase
		para.PC =func_target(time, para);

		// get data for f adjustment
		for(int i= 0; i < para.Ncx; i++) RATES_Fadj.at(i, t)= para.nuCx[i];
		for(int i= 0; i < para.Ngc; i++) INPUTS_Fadj_GC.at(i, t)= para.hGC[i];
		tt++;
	}
	adjust_thresholds_gains(INPUTS_Fadj_GC, RATES_Fadj, para);
	alg_steps_file << "\t \t \t End of coding adjustment phase, T =  "<< para.Tmax_Fadj*para.dt  << "\t Ttot =  "<< time << std::endl;

	//	*************************
	// 	second transient phase (always open loop / driven)
	//	*************************
	// set loop parameter
	para.closedloop=0;		// open loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t Second transient phase"<< "\t RKflag = "<< para.RKflag << "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	// if population activity from this phase is needed
	if(para.intrfflag==2){
		for(int t= 0; t < para.Tmax_trans; t++){

			if(para.taugcmodeflag==0){
				status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, tt*para.dt, hcx);
				update_rates(time, hcx, para);
			}
			else if(para.taugcmodeflag>0){
				status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, tt*para.dt, X);
				update_rates(time, X, para);
			}
			if (status != GSL_SUCCESS){
				printf ("Error in second transient phase, return value=%d\n", status);
				break;
			}

			// set activity of output equal to target to ensure continuity with subsequent phase
			para.PC =func_target(time, para);

			if(para.intrfflag==2){
				// get data for FT
				for(int i= 0; i < para.Ngc; i++) RATES_FTtrans.at(i, t)= para.nuGC[i];
			}

			if(para.learnflag<99) write_single(single_ctxunit_files, single_gc_files, single_output_files, time, para);
			tt++;
		}
	}
	else{
		tt=tt+para.Tmax_trans;
		if(para.taugcmodeflag==0){
			status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, (tt-1)*para.dt, hcx);
			update_rates(time, hcx, para);
		}
		else if(para.taugcmodeflag>0){
			status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, (tt-1)*para.dt, X);
			update_rates(time, X, para);
		}
		if (status != GSL_SUCCESS){
			printf ("Error in second transient phase, return value=%d\n", status);
		}
	}
	alg_steps_file << "\t \t \t End of second transient phase, T =  "<< para.Tmax_trans*para.dt << "\t Ttot =  "<< time << std::endl;

	//	*************************
	//	Fourier analysis for generation of GC interference signals
	//	*************************
	if(para.intrfflag==2){
		alg_steps_file << "\t \t \t Pre-Fourier analysis" << std::endl;
		FTcomp FT_TRANS_comp_ctx;
		calc_FT(RATES_FTtrans, FT_TRANS_comp_ctx, para, para.Tmax_trans,"cereb");
		alg_steps_file << "\t \t \t Calculate artificial overtones" << std::endl;
		calc_artificial_overtones(FT_TRANS_comp_ctx, para, np, json_measures);
	}

	//	*************************
	// 	Save network state
	//	*************************
	for(int i= 0; i < para.Ncx; i++) {
		hcx_save[i]=hcx[i];
	}
	for(int i= 0; i < para.Ngc+para.Ncx; i++) {
		X_save[i]=X[i];
	}
	//	*************************
	//	LS collection phase (always open loop / driven)
	//	*************************
	if(para.LRflag==2 || para.LRflag==3 || para.LRflag==4) 	{
		// set loop parameter
		para.closedloop=0;		// open loop
		para.testphaseflag=0;	// no noise
		para.RKflag=1;			// Runke-Kutta integration
		alg_steps_file << "\t \t \t LS training collection phase"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
		for(int t= 0; t < para.Tmax_learn; t++){

			if(para.taugcmodeflag==0){
				status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, tt*para.dt, hcx);
				update_rates(time, hcx, para);
			}
			else if(para.taugcmodeflag>0){
				status = gsl_odeiv2_driver_apply (RNNGCdriver, &time, tt*para.dt, X);
				update_rates(time, X, para);
			}
			if (status != GSL_SUCCESS){
				printf ("Error in LS collection phase, return value=%d\n", status);
				break;
			}

			// set activity of output equal to target to ensure continuity with subsequent phase
			para.PC =func_target(time, para);

			// collect activity for LS learning
			for(int i= 0; i < para.Ngc; i++) RATES_LS_GC.at(i, t)= para.nuGC[i];

			// write stuff to files
			if ((t % para.write_every)==0 && para.learnflag<99) {

				for(int i= 0; i < para.Ncx; i++)	nuCx_learn_file << para.nuCx[i] << " ";
				nuCx_learn_file << std::endl;

				for(int i= 0; i < para.Ngc; i++) 	nuGC_learn_file << para.nuGC[i] << " ";
				nuGC_learn_file << std::endl;
			}

			if(para.learnflag<99) write_single(single_ctxunit_files, single_gc_files, single_output_files, time, para);
			tt++;
		}
		alg_steps_file << "\t \t \t End of LS training collection phase, T =  "<< para.Tmax_learn*para.dt  << "\t Ttot =  "<< time << std::endl;
	//	*************************
	//	LS/RR offline learning
	//	*************************
		alg_steps_file << "\t \t \t Calculate LS solution"<< std::endl;
		arma::uword rnk;
		// LS learning with GC activity
		LS_solution(RATES_LS_GC, para.Wout_PC, RMSE_train, rnk, para, para.FACsigmaLS_gc);
		//~ json_measures["rank_cereb"].push_back(rnk);

		// calculate training RMSE
		if((para.learnflag>0) && (para.learnflag<99)){
			std::cout << "Training error LS cereb = " << RMSE_train << std::endl;
		}
	}

	gsl_odeiv2_driver_free (RNNdriver);
	gsl_odeiv2_driver_free (RNNGCdriver);

//	**************************************************
//	cereb learning
//	**************************************************
	para.algswitch=2;
	alg_steps_file << "\t \t algswitch = " << para.algswitch<< "\t cerebellar learning"<< std::endl;

	if(para.taugcmodeflag>0){
		for(int i= 0; i < para.Ngc; i++) {
			para.hGC[i]=X_save[i+para.Ncx];
		}
	}

//	*************************
//	RLS/Delta online learning phase (always closed loop)
//	*************************
	if(para.LRflag<2){
		// set loop parameter
		para.closedloop=1;		// closed loop
		para.testphaseflag=0;	// no noise
		if(para.LRflag=1) 	para.RKflag=1;		// Runke-Kutta integration for DELTA
		else 				para.RKflag=0;		// Euler integration for RLS
		alg_steps_file << "\t \t \t RLS/DELTA learning phase"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
		RMSE_train=para.error_signed_train=0;
		stats_rmse.reset();

		// initialise weights
		para.Wout_PC= arma::zeros<arma::vec>(para.Ngc);
		para.Wout_PC.randn();
		//~ para.Wout_PC=para.Wout_PC/sqrt(para.Ngc);

		gsl_odeiv2_driver * RNNGCdriver_learn;

		// DELTA learning
		if(para.LRflag==1) {
			// set initial learning rates for DELTA learning

			std::vector<double> hcxW(para.Ncx+para.Ngc);			// for DELTA learning with instantaneous GCs
			std::vector<double> XW(para.Ncx+para.Ngc+para.Ngc);		// for DELTA learning with integrated GCs

			tt=tt+para.Tmax_learn;
			// DELTA learning for instantaneous GCs
			if(para.taugcmodeflag==0){
				size_t nrnn_delta = static_cast<size_t>(para.Ncx+para.Ngc);
				gsl_odeiv2_system RNNGCdeltasys = {RK_cortex_cereb_weights_rhs, NULL, nrnn_delta, &para};
				RNNGCdriver_learn = gsl_odeiv2_driver_alloc_y_new (&RNNGCdeltasys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

				// initial state
				for(int i= 0; i < para.Ncx; i++) hcxW[i]=hcx[i];
				for(int i= 0; i < para.Ngc; i++) hcxW[i+para.Ncx]=para.Wout_PC[i];	// initial readout weights

				status = gsl_odeiv2_driver_apply (RNNGCdriver_learn, &time, (tt-1)*para.dt, hcxW.data());
				update_rates(time, hcxW.data(), para);
				if (status != GSL_SUCCESS){
					printf ("Error during learning phase, return value=%d\n", status);
					break;
				}
				for(int i= 0; i < para.Ncx; i++) hcx[i]=hcxW[i];
				for(int i= 0; i < para.Ngc; i++) para.Wout_PC[i]=hcxW[i+para.Ncx];
			}
			else if(para.taugcmodeflag>0){
				size_t nrnn_delta = static_cast<size_t>(para.Ncx+2*para.Ngc);
				gsl_odeiv2_system RNNGCdeltasys = {RK_cortex_cereb_dyn_weights_rhs, NULL, nrnn_delta, &para};
				RNNGCdriver_learn = gsl_odeiv2_driver_alloc_y_new (&RNNGCdeltasys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

				// initial state
				for(int i= 0; i < para.Ncx+para.Ngc; i++)	XW[i]=X[i];
				for(int i= 0; i < para.Ngc; i++) 			XW[i+para.Ncx+para.Ngc]=para.Wout_PC[i];	// initial readout weights

				status = gsl_odeiv2_driver_apply (RNNGCdriver_learn, &time, (tt-1)*para.dt, XW.data());
				update_rates(time, XW.data(), para);
				if (status != GSL_SUCCESS){
					printf ("Error during learning phase, return value=%d\n", status);
					break;
				}
				for(int i= 0; i < para.Ncx+para.Ngc; i++)	X[i]=XW[i];
				for(int i= 0; i < para.Ncx; i++) 			hcx[i]=XW[i];
				for(int i= 0; i < para.Ngc; i++) 			para.hGC[i]=XW[i+para.Ncx];
				for(int i= 0; i < para.Ngc; i++) 			para.Wout_PC[i]=XW[i+para.Ncx+para.Ngc];
			}

			RMSE_train=NAN;
			alg_steps_file << "\t \t \t End of DELTA learning phase, T =  "<< para.Tmax_learn*para.dt << "\t Ttot =  "<< time << std::endl;

			gsl_odeiv2_driver_free (RNNGCdriver_learn);
		}
		// RLS learning
		else{
			// set up initial inverse correlation matrix
			RLS_cereb.P=arma::zeros<arma::mat>(para.Ngc,para.Ngc);
			for(int i= 0; i < para.Ngc; i++) {
				RLS_cereb.P(i,i)=1.0/para.alpha;
			}

			// 	set time constants
			set_taus(para,para.dt);
			// compute output
			para.PC=output_rate(para.Wout_PC, para.nuGC);
			for(int t= 0; t < para.Tmax_learn; t++){
				time=tt*para.dt;
				// instantaneous GCs
				if(para.taugcmodeflag==0){
					euler_step_cortex_cerebellum(time, hcx, para);
				}
				// integrated GCs
				else if(para.taugcmodeflag>0){
					euler_step_cortex_cerebellum_dyn(time, X, para);
				}

				if ((t % para.RLSlearn_every)==0) {
					RLS_step_cerebellum_allGCs(RLS_cereb, time, para);
				}
				stats_rmse(para.error_signed_train*para.error_signed_train);

				// compute output
				para.PC=output_rate(para.Wout_PC, para.nuGC);
				OUTvec[t]=para.PC;
				Target[t]=func_target(time, para);

				// write stuff to files
				if ((t % para.write_every)==0 && para.learnflag<99) {
					for(int i= 0; i < para.Ncx; i++)	nuCx_learn_file << para.nuCx[i] << " ";
					nuCx_learn_file << std::endl;

					for(int i= 0; i < para.Ngc; i++)	nuGC_learn_file << para.nuGC[i] << " ";
					nuGC_learn_file << std::endl;
				}

				if (para.learnflag<99) {
					if ((t % para.write_every)==0) {
						LR_file << RLS_cereb.LR << std::endl;
						error_signed_train_file << para.error_signed_train << std::endl;
						error_full_train_file << sqrt(stats_rmse.mean()) << std::endl;
					}
				}
				if(para.learnflag<99) write_single(single_ctxunit_files, single_gc_files, single_output_files, time, para);
				tt++;
			}
			alg_steps_file << "\t \t \t End of RLS learning phase, T =  "<< para.Tmax_learn*para.dt << "\t Ttot =  "<< time << std::endl;

			// calculate training RMSE
			if(para.learnflag>0){
				RMSE_train=sqrt(arma::mean((OUTvec-Target)%(OUTvec-Target)));
				if(para.learnflag<99)	std::cout << "Training error RLS/DELTA = " << RMSE_train << std::endl;
			}
		}

	}
//	*************************
//	test phase (always closed loop)
//	*************************

	gsl_odeiv2_driver * RNNGCdriver_test = gsl_odeiv2_driver_alloc_y_new (&RNNGCsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

	// set loop parameter
	para.closedloop=para.loop_test_flag; 	// closed/open loop set by loop_test_flag
	para.testphaseflag=1;					// testphase: add noise
	if(para.LRflag<2) 	para.RKflag=0;		// Euler integration for RLS variants
	else 				para.RKflag=1;		// Runke-Kutta integration for LS variants
	alg_steps_file << "\t \t \t Test phase"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	RMSE_test=0;
	int ti=0;
	// Euler integration (with noise)
	if(para.RKflag==0){
		// 	set time constants
		set_taus(para,para.dt_euler);
		for(int t= 0; t < para.Tmax_test_euler; t++){
			bool bool_traj=( (t % int(para.dt/para.dt_euler))==0);

			time+=para.dt_euler;
			if(para.taugcmodeflag==0){
				euler_step_cortex_cerebellum(time, hcx, para);
			}
			else if(para.taugcmodeflag>0){
				euler_step_cortex_cerebellum_dyn(time, X, para);
			}
			// compute output
			para.PC=output_rate(para.Wout_PC, para.nuGC);

			// compute test error
			error_signed_test = para.PC-func_target(time, para);
			RMSE_test+=error_signed_test*error_signed_test;

			// get open closed trajectory
			if( (ti < para.Tmax_FM) && bool_traj) {
				if(para.taugcmodeflag==0){
					for(int i= 0; i < para.Ncx; i++) 			INPUTS_CLOSED.at(i, ti)= hcx[i];
				}
				else if(para.taugcmodeflag>0){
					for(int i= 0; i < para.Ncx+para.Ngc; i++) 	INPUTS_CLOSED.at(i, ti)= X[i];
				}
				ti++;
			}

			// write stuff to files
			if ((t % (para.write_every* int(para.dt/para.dt_euler)))==0) {
				if(para.learnflag<99){
					error_signed_test_file << error_signed_test << std::endl;

					for(int i= 0; i < para.Ncx; i++)	nuCx_test_file << para.nuCx[i] << " ";
					nuCx_test_file << std::endl;
					for(int i= 0; i < para.Ngc; i++) 	nuGC_test_file << para.nuGC[i] << " ";
					nuGC_test_file << std::endl;
				}
				if(para.ROwriteflag==1) {
					PC_test_files[np] << para.PC << "\t";
				}
			}

			if(para.learnflag<99) write_single(single_ctxunit_files, single_gc_files, single_output_files, time, para);
			tt++;
		}
		// calculate test RMSE
		RMSE_test=sqrt(RMSE_test/double(para.Tmax_test_euler))/sqrt(para.signalvar);
		alg_steps_file << "\t \t \t End of test phase, T =  "<< para.Tmax_test_euler*para.dt_euler << "\t Ttot =  "<< time << std::endl;
	}
	// RK integration (noise-free)
	else{
		for(int t= 0; t < para.Tmax_test; t++){

			if(para.taugcmodeflag==0){
				status = gsl_odeiv2_driver_apply (RNNGCdriver_test, &time, tt*para.dt, hcx);
				update_rates(time, hcx, para);
			}
			else if(para.taugcmodeflag>0){
				status = gsl_odeiv2_driver_apply (RNNGCdriver_test, &time, tt*para.dt, X);
				update_rates(time, X, para);
			}
			if (status != GSL_SUCCESS){
				printf ("Error in test phase, return value=%d\n", status);
				break;
			}

			// compute output
			para.PC=output_rate(para.Wout_PC, para.nuGC);

			// compute test error
			error_signed_test = para.PC-func_target(time, para);
			RMSE_test+=error_signed_test*error_signed_test;

			// get closed loop trajectory
			if(t < para.Tmax_FM) {
				if(para.taugcmodeflag==0){
					for(int i= 0; i < para.Ncx; i++) INPUTS_CLOSED.at(i, t)= hcx[i];
				}
				else if(para.taugcmodeflag>0){
					for(int i= 0; i < para.Ncx+para.Ngc; i++) INPUTS_CLOSED.at(i, t)= X[i];
				}
			}

			// write stuff to files
			if ((t % para.write_every)==0) {
				if(para.learnflag<99){
					error_signed_test_file << error_signed_test << std::endl;

					for(int i= 0; i < para.Ncx; i++)	nuCx_test_file << para.nuCx[i] << " ";
					nuCx_test_file << std::endl;
					for(int i= 0; i < para.Ngc; i++) 	nuGC_test_file << para.nuGC[i] << " ";
					nuGC_test_file << std::endl;
				}
				if(para.ROwriteflag==1) {
					PC_test_files[np] << para.PC << "\t";
				}
			}

			if(para.learnflag<99) write_single(single_ctxunit_files, single_gc_files, single_output_files, time, para);
			tt++;
		}
		// calculate test RMSE
		RMSE_test=sqrt(RMSE_test/double(para.Tmax_test))/sqrt(para.signalvar);
		alg_steps_file << "\t \t \t End of test phase, T =  "<< para.Tmax_test*para.dt << "\t Ttot =  "<< time << std::endl;
	}

	if(para.learnflag<99) {
		std::cout << "Test error = " << RMSE_test << std::endl;
	}
	if(para.ROwriteflag==1) {
	 	PC_test_files[np] << std::endl;
	}

	gsl_odeiv2_driver_free (RNNGCdriver_test);

//	**************************************************
// 	Calculate Floquet multipliers
//	**************************************************
	calc_floquet_RK_vector(X_save, para, json_measures, alg_steps_file);

//	*************************
// 	Calculate norm of output weights
//	*************************
	double Wnorm_cereb=calc_wnorm(para.Wout_PC);

//	**************************************************
// 	Calculate Fourier components
//	**************************************************
	gsl_odeiv2_driver * RNNGCdriver_FT = gsl_odeiv2_driver_alloc_y_new (&RNNGCsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

	// set loop parameter
	para.closedloop=0;		// open loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t Transient phase for FT calculations"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	// set initial network state
	for(int i= 0; i < para.Ncx; i++) {
		hcx[i]=hcx_save[i];
	}
	for(int i= 0; i < para.Ncx+para.Ngc; i++) {
		X[i]=X_save[i];
	}
	stats_hcx.reset(); 	stats_nuCx.reset();
	stats_hgc.reset(); 	stats_GC.reset();
	time=0;

	tt=para.Tmax_trans;
	if(para.taugcmodeflag==0){
		status = gsl_odeiv2_driver_apply (RNNGCdriver_FT, &time, (para.Tmax_trans-1)*para.dt, hcx);
		update_rates(time, hcx, para);
	}
	else if(para.taugcmodeflag>0){
		status = gsl_odeiv2_driver_apply (RNNGCdriver_FT, &time, (para.Tmax_trans-1)*para.dt, X);
		update_rates(time, X, para);
	}
	if (status != GSL_SUCCESS){
		printf ("Error during Fourier transient phase, return value=%d\n", status);
		break;
	}
	alg_steps_file << "\t \t \t End of FT transient phase, T =  "<< para.Tmax_trans*para.dt << "\t Ttot =  "<< time << std::endl;

	alg_steps_file << "\t \t \t FT trajectory extraction phase"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	// multiple periods for FT trajectory
	for(int t= 0; t < para.Tmax_FT; t++){
		fsparse_pop[t]=0;
		nuGC_avrg_pop[t]=0;

		if(para.taugcmodeflag==0){
			status = gsl_odeiv2_driver_apply (RNNGCdriver_FT, &time, tt*para.dt, hcx);
			update_rates(time, hcx, para);
		}
		else if(para.taugcmodeflag>0){
			status = gsl_odeiv2_driver_apply (RNNGCdriver_FT, &time, tt*para.dt, X);
			update_rates(time, X, para);
			for(int i= 0; i < para.Ncx; i++) {
				hcx[i]=X[i];
			}
		}
		if (status != GSL_SUCCESS){
			printf ("Error during FT phase, return value=%d\n", status);
			break;
		}

		// get data for FT
		for(int i= 0; i < para.Ncx; i++) RATES_FT.at(i, t)= para.nuCx[i];
		for(int i= 0; i < para.Ngc; i++) RATES_FT_GC.at(i, t)= para.nuGC[i];

		// get open loop trajectory
		if(t < para.Tmax_FM) {
			if(para.taugcmodeflag==0){
				for(int i= 0; i < para.Ncx; i++) INPUTS_OPEN.at(i, t)= hcx[i];
			}
			else if(para.taugcmodeflag>0){
				for(int i= 0; i < para.Ncx+para.Ngc; i++) INPUTS_OPEN.at(i, t)= X[i];
			}
		}

		// set activity of output equal to target
		para.PC =func_target(time, para);

		for(int i= 0; i < para.Ncx; i++){
			stats_hcx(hcx[i]);
			stats_nuCx(para.nuCx[i]);
		}
		// cerebellar population
		for(int i= 0; i < para.Ngc; i++){
			stats_hgc(para.hGC[i]);	stats_GC(para.nuGC[i]);
			if(std::fabs(para.nuGC[i])>0.01)	{
				fsparse[i]+=1.0;
				fsparse_pop[t]+=1.0;
			}
			nuGC_avrg_pop[t]+=para.nuGC[i];
		}

		// write to files
		if(para.learnflag<99){
			if ((t % para.write_every)==0) {
				for(int i= 0; i < para.Ncx; i++)	nuCx_FT_file << para.nuCx[i] << " ";
				nuCx_FT_file << std::endl;

				for(int i= 0; i < para.Ngc; i++) 	nuGC_FT_file << para.nuGC[i] << " ";
				nuGC_FT_file << std::endl;

				for(int i= 0; i < para.Ngc; i++) 	hGC_FT_file << para.hGC[i] << " ";
				hGC_FT_file << std::endl;
			}
		}

		if(para.learnflag<99) write_single(single_ctxunit_FT_files, single_gc_FT_files, single_output_FT_files, time, para);
		tt++;
	}
	alg_steps_file << "\t \t \t End of FT trajectory extraction phase, T =  "<< para.Tmax_FT*para.dt  << "\t Ttot =  "<< time << std::endl;

	gsl_odeiv2_driver_free (RNNGCdriver_FT);

//	*************************
// 	Calculate measures of interest
//	*************************
	// calculate life-time sparsity stats
	fsparse=fsparse/double(para.Tmax_FT);
	double fmean=arma::mean(fsparse);
	double fstd=arma::stddev(fsparse);
	double fmin=arma::min(fsparse);
	double fmax=arma::max(fsparse);

	// calculate population sparsity stats
	fsparse_pop=fsparse_pop/double(para.Ngc);
	double fpop_mean=arma::mean(fsparse_pop);
	double fpop_std=arma::stddev(fsparse_pop);
	double fpop_min=arma::min(fsparse_pop);
	double fpop_max=arma::max(fsparse_pop);

	// calculate fluctuation of population mean
	nuGC_avrg_pop=nuGC_avrg_pop/double(para.Ngc);
	double nuGC_avrg_std=arma::stddev(nuGC_avrg_pop);
	double nuGC_avrg_cv=nuGC_avrg_std/arma::mean(nuGC_avrg_pop);

	// write stats to file
	std::ofstream stats_file(para.folder+"stats.dat");		if (!stats_file) exit(EXIT_FAILURE);

	json_stats["hcx_mean"].push_back(stats_hcx.mean());
	json_stats["hcx_std"].push_back(stats_hcx.stddev());
	json_stats["nucx_mean"].push_back(stats_nuCx.mean());
	json_stats["nucx_std"].push_back(stats_nuCx.stddev());
	json_stats["hgc_mean"].push_back(stats_hgc.mean());
	json_stats["hgc_std"].push_back(stats_hgc.stddev());
	json_stats["gc_mean"].push_back(stats_GC.mean());
	json_stats["gc_std"].push_back(stats_GC.stddev());
	json_stats["gc_pop_std"].push_back(nuGC_avrg_std);
	json_stats["gc_pop_cv"].push_back(nuGC_avrg_cv);
	json_stats["f_mean"].push_back(fmean);
	json_stats["f_std"].push_back(fstd);
	json_stats["f_min"].push_back(fmin);
	json_stats["f_max"].push_back(fmax);
	json_stats["fpop_mean"].push_back(fpop_mean);
	json_stats["fpop_std"].push_back(fpop_std);
	json_stats["fpop_min"].push_back(fpop_min);
	json_stats["fpop_max"].push_back(fpop_max);
	json_stats["thres_mean"].push_back(arma::mean(para.threshold_vector));
	json_stats["gain_mean"].push_back(arma::mean(para.gain_vector));

	stats_file << std::setw(4) << json_stats << '\n';
	stats_file.close();

	// calculate Fourier components & FT dimensionality
	alg_steps_file << "\t \t \t Fourier analysis" << std::endl;
	FTcomp FTcomp_ctx;
	FTcomp FTcomp_cereb;
	calc_FT(RATES_FT, FTcomp_ctx, para, para.Tmax_FT,"ctx");
	calc_FT(RATES_FT_GC, FTcomp_cereb, para, para.Tmax_FT,"cereb");

	// calculate dimensionality
	alg_steps_file << "\t \t \t Calculate dimensionality" << std::endl;
	double DIM_ctx=calc_dim_noEV(RATES_FT,para.Tmax_FT);
	double DIM_cereb=calc_dim_noEV(RATES_FT_GC,para.Tmax_FT);

	// write to files: ctx
	alg_steps_file << "\t \t \t Write population measures to files"<< std::endl;
	json_measures["DIM_ctx"].push_back(DIM_ctx);
	json_measures["DIM_FT_ctx"].push_back(FTcomp_ctx.dim_signal);
	json_measures["COND_ctx"].push_back(FTcomp_ctx.cond_signal);
	json_measures["COND_total_ctx"].push_back(FTcomp_ctx.cond_total);

	// write to files: cereb
	json_measures["DIM_cereb"].push_back(DIM_cereb);
	json_measures["DIM_FT_cereb"].push_back(FTcomp_cereb.dim_signal);
	json_measures["COND_cereb"].push_back(FTcomp_cereb.cond_signal);
	json_measures["COND_total_cereb"].push_back(FTcomp_cereb.cond_total);

	// RMSE of driven and closed loop trajectories
	arma::vec DEVT(para.Tmax_FM);
	arma::vec OPENR(para.Tmax_FM);
	for(int t= 0; t < para.Tmax_FM; t++){
		DEVT(t)=arma::norm(INPUTS_OPEN.col(t)-INPUTS_CLOSED.col(t));
		OPENR(t)=arma::norm(INPUTS_OPEN.col(t));
	}
	//~ double rmse_traj= arma::mean(DEVT)/sqrt(INPUTS_OPEN.n_rows);
	double rmse_traj= arma::mean(DEVT)/arma::mean(OPENR);
	json_measures["Traj_rmse_"].push_back(rmse_traj);

//	*************************
// 	write output to files
//	*************************
	para.folder=para.path+"parameter_"+ std::to_string (np+1) + "/";

	if(para.learnflag<99){
		alg_steps_file << "\t \t \t Write single trial data to files"<< std::endl;
		//	connectivity files
		std::ofstream Wcx_file(para.path_connect+"Wcx.dat");				if (!Wcx_file) exit(EXIT_FAILURE);
		std::ofstream Wout_file(para.path_connect+"Wout.dat");				if (!Wout_file) exit(EXIT_FAILURE);
		std::ofstream Wfb_file(para.path_connect+"Wfb.dat");				if (!Wfb_file) exit(EXIT_FAILURE);
		std::ofstream Wfb_interf_file(para.path_connect+"Wfb_interf.dat");	if (!Wfb_interf_file) exit(EXIT_FAILURE);
		std::ofstream Wmf_file(para.path_connect+"Wmf.dat");				if (!Wmf_file) exit(EXIT_FAILURE);
		std::ofstream Wout_PC_file(para.path_connect+"J.dat");					if (!Wout_PC_file) exit(EXIT_FAILURE);
		std::ofstream MFidx_file(para.path_connect+"MFidx.dat");			if (!MFidx_file) exit(EXIT_FAILURE);

		//	write stuff to connectivity-files
		for(int j= 0; j < para.Ngc; j++){
			Wout_PC_file << para.Wout_PC[j] << "\t";
		}
		Wout_PC_file << std::endl;
		for(int j= 0; j < para.Ncx; j++) {
			Wout_file << para.Wout[j] << "\t";
		}
		Wout_file << std::endl;
		for(int i= 0; i < para.Ncx; i++) {
			for(int j= 0; j < para.Wcx[i].Nsyn; j++) 	Wcx_file  << i << "\t" << para.Wcx[i].idx[j] << "\t" << para.Wcx[i].weight[j] << std::endl;
			Wfb_file  << para.Wfb[i] << std::endl;
			if(para.intrfflag==1){
				for(int k=0; k < para.Ncomp_intf; k++) 	Wfb_interf_file  << para.Wfb_interf.at(i,k) << "\t";
				Wfb_interf_file << std::endl;
			}
			else if(para.intrfflag==2){
				for(int j= 0; j < para.Ncomp_intf; j++) {
					Wfb_interf_file  << para.Intf_weights(i,j) << "\t";
				}
				Wfb_interf_file  <<  std::endl;
			}
		}
		for(int i= 0; i < para.Ngc; i++) {
			for(int j= 0; j < para.Nsyn; j++) {
				Wmf_file << i+1 << "\t" << para.Wmf[i].idx[j]+1 << "\t" << para.Wmf[i].weight[j]<< std::endl;
			}
		}

		for(int i= 0; i < para.Ncx; i++) {
			json_trial["tau_distrb"].push_back(para.tau_distrb[i]);
		}
		for(int i= 0; i < para.Ngc; i++) {
			json_trial["theta"].push_back(para.threshold_vector[i]);
			json_trial["gain"].push_back(para.gain_vector[i]);
			json_trial["tauGC_distrb"].push_back(para.tauGC_distrb[i]);
			json_trial["Wgoc_distrb"].push_back(para.Wgoc_distrb[i]);
		}

		for (int j= 0;j<para.Nmf;j++) MFidx_file << j+1 << "\t" << para.Ncx_to_MF_idx[j]+1 << "\t" << 1 << std::endl;

		FTcomp_ctx.Angles_all.save(angles_ctx_file,arma::raw_ascii);
		FTcomp_cereb.Angles_all.save(angles_cereb_file,arma::raw_ascii);

		FTcomp_ctx.S_total.save(S_total_ctx_file,arma::raw_ascii);
		FTcomp_cereb.S_total.save(S_total_cereb_file,arma::raw_ascii);

		// write to file
		std::ofstream single_trial_file(para.folder+"single_trial.dat");		if (!single_trial_file) exit(EXIT_FAILURE);
		single_trial_file << std::setw(4) << json_trial << '\n';
		single_trial_file.close();

		Wcx_file.close();
		Wmf_file.close();
		Wout_file.close();
		Wfb_file.close();
		Wfb_interf_file.close();
		Wout_PC_file.close();
		MFidx_file.close();
	}

	bool write_FULL=1;
	alg_steps_file << "\t \t \t Write cereb summary data to files"<< std::endl;
	write_measures_final(json_measures, FTcomp_cereb, RMSE_train, RMSE_test, Wnorm_cereb, "cereb", para, write_FULL);

	simcount_file << 1 << std::endl;

	} // end of Nreal-loop

	PC_test_files[np].close();

	} // end of Npara-loop

//	*************************
// 	terminate program
	alg_steps_file << "Terminate"<< std::endl;
   	nuCx_learn_file.close();
   	nuCx_test_file.close();
	nuCx_FT_file.close();

	nuGC_learn_file.close();
	nuGC_test_file.close();
	nuGC_FT_file.close();
	hGC_FT_file.close();
	Inputs_file.close();

	error_signed_test_file.close();
	error_signed_train_file.close();
	error_full_train_file.close();
	LR_file.close();

	single_ctxunit_files.close();
	single_gc_files.close();
	single_output_files.close();
	single_ctxunit_FT_files.close();
	single_gc_FT_files.close();
	single_output_FT_files.close();

	freq_file.close();
	amp_file.close();
	phase_file.close();
	freq_intf_file.close();
	amp_intf_file.close();
	phase_intf_file.close();
	angles_ctx_file.close();
	angles_cereb_file.close();
	S_total_ctx_file.close();
	S_total_cereb_file.close();

	alg_steps_file.close();
	simcount_file.close();

	gsl_rng_free (para.r);

  return(0);

}
