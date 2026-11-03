// added boost::program_options

#include <vector>
#include <armadillo>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <sys/stat.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "headers/json.hpp"
#include "headers/functions.h"
#include "headers/allstructs_v5.h"

#include "headers/filter_GCs_v2.h"
#include "headers/target_signal_v6.h"
#include "headers/euler_updates_v1.h"
#include "headers/RK_rhs_v1.h"
#include "headers/RK_weights_rhs_v1.h"
#include "headers/update_rates_ctx_v1.h"
#include "headers/RLS_step_cortex_v3.h"
#include "headers/calc_FT_v7.h"
#include "headers/calc_floquet_RK_vector_v1.h"
#include "headers/calc_connectivityEVs_v2.h"
#include "headers/LS_solution_v4.h"
#include "headers/generate_weights_v3.h"
#include "headers/output_calculations_v2.h"
#include "headers/set_scan_variables.h"
#include "headers/set_parameters_flags_ctx_only_v2.h"
#include "headers/write_parameters_v2.h"
#include "headers/write_measures_final_v4.h"
#include "headers/set_taus_v1.h"
#include "headers/write_single_ctx_v1.h"

int main (int argc, char* argv[]) {

//------------------------------------------------------------------------
//	set parameters, flags, write to file & print some output
//------------------------------------------------------------------------
	struct netpara para;
	para.version="cortex_only_v06";

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

	std::ofstream Z_file(para.path+"Z.dat");						if (!Z_file) exit(EXIT_FAILURE);
	std::ofstream Inputs_file(para.path+"Inputs.dat");				if (!Inputs_file) exit(EXIT_FAILURE);

	std::ofstream freq_file(para.path+"frequency.dat");			if (!freq_file) exit(EXIT_FAILURE);
	std::ofstream amp_file(para.path+"amp.dat");				if (!amp_file) exit(EXIT_FAILURE);
	std::ofstream phase_file(para.path+"phase.dat");			if (!phase_file) exit(EXIT_FAILURE);
	std::ofstream freq_intf_file(para.path+"frequency_intf.dat");		if (!freq_intf_file) exit(EXIT_FAILURE);
	std::ofstream amp_intf_file(para.path+"amp_intf.dat");			if (!amp_intf_file) exit(EXIT_FAILURE);
	std::ofstream phase_intf_file(para.path+"phase_intf.dat");		if (!phase_intf_file) exit(EXIT_FAILURE);

	std::ofstream angles_ctx_file(para.path+"angles_ctx.dat");		if (!angles_ctx_file) exit(EXIT_FAILURE);
	std::ofstream S_total_ctx_file(para.path+"stotal_ctx.dat");		if (!S_total_ctx_file) exit(EXIT_FAILURE);

	std::ofstream error_signed_train_file(para.path+"error_signed_train.dat");	if (!error_signed_train_file) exit(EXIT_FAILURE);
	std::ofstream error_full_train_file(para.path+"error_full_train.dat");		if (!error_full_train_file) exit(EXIT_FAILURE);
	std::ofstream error_signed_test_file(para.path+"error_test.dat");			if (!error_signed_test_file) exit(EXIT_FAILURE);
	std::ofstream LR_file(para.path+"LR.dat");									if (!LR_file) exit(EXIT_FAILURE);

	std::ofstream simcount_file(para.path+"simcount.dat");				if (!simcount_file) exit(EXIT_FAILURE);

	std::ofstream single_ctxunit_files(para.path+"single_ctxunit.dat");	if (!single_ctxunit_files) exit(EXIT_FAILURE);
	std::ofstream single_output_files(para.path+"single_output.dat");	if (!single_output_files) exit(EXIT_FAILURE);
	std::ofstream single_ctxunit_FT_files(para.path+"single_ctxunit_FT.dat");	if (!single_ctxunit_FT_files) exit(EXIT_FAILURE);
	std::ofstream single_output_FT_files(para.path+"single_output_FT.dat");	if (!single_output_FT_files) exit(EXIT_FAILURE);

	std::ofstream alg_steps_file(para.path+"alg_steps.dat");			if (!alg_steps_file) exit(EXIT_FAILURE);

	std::ofstream Z_test_files[para.Npara];

	for(int np= 0; np < para.Npara; np++) {
		para.folder=para.path+"parameter_"+ std::to_string (np+1) + "/";

		Z_test_files[np].open(para.folder+"Z_test.dat", std::ofstream::out | std::ofstream::trunc);
		Z_test_files[np].close();
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

		Z_test_files[np].open(para.folder+"Z_test.dat",std::ios::app);

		// set parameters
		if(para.learnflag==100){
			para.sigmaN=SCANPARA(np);
			para.sigmaNgc=para.sigmaN;
			std::cout << "sigmaN = " << para.sigmaN << std::endl;
		}
		else if(para.learnflag==102){
			para.Ncx=SCANPARA(np);
			para.Nmf=para.Ncx;
			std::cout << "Ncx = " << para.Ncx << std::endl;
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
			std::cout << "eta_ini = " << para.delta_eta << std::endl;
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
		else if(para.learnflag==128){
			para.rho_Wfb12=SCANPARA(np);
			para.rho_Wfb={ 	{1, 				1, 	para.rho_Wfb12}, 				// 3x3
							{1, 				1, 	1},
							{para.rho_Wfb12, 	1, 	1} };
			std::cout << "rho_Wfb_12 = " << para.rho_Wfb12 << std::endl;
		}
		else if(para.learnflag==129){
			para.sigAmp_intf=SCANPARA(np);
			std::cout << "sigAmp_intf = " << para.sigAmp_intf << std::endl;
		}

	para.dt=para.DT/para.lambda;
	para.dt_euler=para.DT_euler/para.lambda;
	para.Tmax_learn = round(para.duration_learn/para.DT);
	para.Tmax_test = round(para.duration_test/para.DT);
	para.Tmax_test_euler = round(para.duration_test/para.DT_euler);
	para.Tmax_FT = round(para.duration_FT/para.DT);
	para.Tmax_trans = round(para.duration_trans/para.DT);
	para.Tmax_FM = round(1.0/para.DT);						// one period for FM input arrays

	int status;

	double hcx[para.Ncx],hcx_save[para.Ncx],hcxini[para.Ncx];
	double time,error_signed_test,RMSE_train,RMSE_test;

	struct RLS RLS_ctx;
	RLS_ctx.N=para.Ncx;
	RLS_ctx.dwout=arma::zeros<arma::vec>(para.Ncx);
	RLS_ctx.P=arma::zeros<arma::mat>(para.Ncx,para.Ncx);

	arma::mat RATES_FT(para.Ncx, para.Tmax_FT);
	arma::mat RATES_LS(para.Ncx, para.Tmax_learn);
	arma::mat INPUTS_CLOSED(para.Ncx, para.Tmax_FM);
	arma::mat INPUTS_OPEN(para.Ncx, para.Tmax_FM);

	arma::colvec nuCx_avrg_pop(para.Tmax_FT);
	arma::colvec OUTvec(para.Tmax_learn);
	arma::colvec Target(para.Tmax_learn);

	arma::running_stat<double> stats_hcx;
	arma::running_stat<double> stats_nuCx;
	arma::running_stat<double> stats_rmse;
//------------------------------------------------------------------------
// 	generate target signal(s)
//------------------------------------------------------------------------
	alg_steps_file << "\t \t write target signal(s) to: "<< para.path+"parameter_"+ std::to_string (np+1) + "/" << std::endl;
	target_signal(para,freq_file,amp_file,phase_file,freq_intf_file,amp_intf_file,phase_intf_file,alg_steps_file);

//	*************************
// 	Set sigmaLS multiplier
//	*************************
	para.FACsigmaLS_ctx=para.sigAmp;

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
  	gsl_odeiv2_system RNNsys = {RK_cortex_rhs, NULL, nrnn, &para};

	//~ gsl_odeiv2_driver * RNNdriver = gsl_odeiv2_driver_alloc_y_new (&RNNsys, gsl_odeiv2_step_rk8pd,para.dt, para.INT_delta_abs,para.INT_delta_rel);
	gsl_odeiv2_driver * RNNdriver = gsl_odeiv2_driver_alloc_y_new (&RNNsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

	//	*************************
	// 	set initial conditions
	//	*************************
	para.nuCx= arma::zeros<arma::vec>(para.Ncx);
	for(int i= 0; i < para.Ncx; i++) {
		hcx[i]=hcxini[i];
	}
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

	//	*************************
	// 	Save network state
	//	*************************
	for(int i= 0; i < para.Ncx; i++) {
		hcx_save[i]=hcx[i];
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

			status = gsl_odeiv2_driver_apply (RNNdriver, &time, tt*para.dt, hcx);
			update_rates_ctx(time, hcx, para);
			if (status != GSL_SUCCESS){
				printf ("Error in LS collection phase, return value=%d\n", status);
				break;
			}

			// set activity of output equal to target to ensure continuity with subsequent phase
			para.Z =func_target(time, para);

			// collect activity for LS learning
			for(int i= 0; i < para.Ncx; i++) RATES_LS.at(i, t)= para.nuCx[i];

			// write stuff to files
			if ((t % para.write_every)==0 && para.learnflag<99) {

				for(int i= 0; i < para.Ncx; i++)	nuCx_learn_file << para.nuCx[i] << " ";
				nuCx_learn_file << std::endl;

				//~ Z_file << para.Z << std::endl;
			}

			if(para.learnflag<99) write_single_ctx(single_ctxunit_files, single_output_files, time, para);
			tt++;
		}
		alg_steps_file << "\t \t \t End of LS training collection phase, T =  "<< para.Tmax_learn*para.dt  << "\t Ttot =  "<< time << std::endl;
	//	*************************
	//	LS/RR offline learning
	//	*************************
		alg_steps_file << "\t \t \t Calculate LS solution"<< std::endl;
		arma::uword rnk;
		// LS learning with cortical activity
		LS_solution(RATES_LS, para.Wout, RMSE_train, rnk, para, para.FACsigmaLS_ctx);
		//~ json_measures["rank_ctx"].push_back(rnk);

		// calculate training RMSE
		if((para.learnflag>0) && (para.learnflag<99)){
			std::cout << "Training error LS ctx = " << RMSE_train << std::endl;
		}
	}

	gsl_odeiv2_driver_free (RNNdriver);

//	**************************************************
// 	Calculate effective connectivity EV spectrum (linear RNN and LS learning only)
//	**************************************************
	calc_connectivityEVs(para, json_measures, alg_steps_file);

//	**************************************************
//	ctx learning
//	**************************************************
	para.algswitch=1;
	alg_steps_file << "\t \t algswitch = " << para.algswitch<< "\t cortical learning"<< std::endl;

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
		para.Wout= arma::zeros<arma::vec>(para.Ncx);
		para.Wout.randn();

		gsl_odeiv2_driver * RNNdriver_learn;

		// DELTA learning
		if(para.LRflag==1) {
			// set initial learning rates for DELTA learning

			std::vector<double> hcxW(para.Ncx+para.Ncx);

			tt=tt+para.Tmax_learn;

			size_t nrnn_delta = static_cast<size_t>(para.Ncx+para.Ncx);
			gsl_odeiv2_system RNNdeltasys = {RK_cortex_weights_rhs, NULL, nrnn_delta, &para};
			RNNdriver_learn = gsl_odeiv2_driver_alloc_y_new (&RNNdeltasys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

			// initial state
			for(int i= 0; i < para.Ncx; i++) hcxW[i]=hcx[i];
			for(int i= 0; i < para.Ncx; i++) hcxW[i+para.Ncx]=para.Wout[i];		// initial readout weights

			status = gsl_odeiv2_driver_apply (RNNdriver_learn, &time, (tt-1)*para.dt, hcxW.data());
			update_rates_ctx(time, hcxW.data(), para);
			if (status != GSL_SUCCESS){
				printf ("Error during learning phase, return value=%d\n", status);
				break;
			}
			for(int i= 0; i < para.Ncx; i++) hcx[i]=hcxW[i];
			for(int i= 0; i < para.Ncx; i++) para.Wout[i]=hcxW[i+para.Ncx];

			RMSE_train=NAN;
			alg_steps_file << "\t \t \t End of DELTA learning phase, T =  "<< para.Tmax_learn*para.dt << "\t Ttot =  "<< time << std::endl;

			gsl_odeiv2_driver_free (RNNdriver_learn);
		}
		// RLS learning
		else{
			// set up initial inverse correlation matrix
			RLS_ctx.P=arma::zeros<arma::mat>(para.Ncx,para.Ncx);
			for(int i= 0; i < para.Ncx; i++) {
				RLS_ctx.P(i,i)=1.0/para.alpha;
			}

			// 	set time constants
			set_taus(para,para.dt);
			// compute output
			para.Z=output_rate(para.Wout, para.nuCx);
			for(int t= 0; t < para.Tmax_learn; t++){
				time=tt*para.dt;
				euler_step_cortex(time, hcx, para);

				if ((t % para.RLSlearn_every)==0) {
					RLS_step_cortex(RLS_ctx, time, para);
				}
				stats_rmse(para.error_signed_train*para.error_signed_train);

				// compute output
				para.Z=output_rate(para.Wout, para.nuCx);
				OUTvec[t]=para.Z;
				Target[t]=func_target(time, para);

				// write stuff to files
				if ((t % para.write_every)==0 && para.learnflag<99) {
					for(int i= 0; i < para.Ncx; i++)	nuCx_learn_file << para.nuCx[i] << " ";
					nuCx_learn_file << std::endl;
					Z_file << para.Z << std::endl;
				}

				if (para.learnflag<99) {
					if ((t % para.write_every)==0) {
						//~ calc_wnorm_change(RLS_ctx, para.Ncx, para);
						LR_file << RLS_ctx.LR << std::endl;
						//~ LR_file << et << std::endl;
						error_signed_train_file << para.error_signed_train << std::endl;
						error_full_train_file << sqrt(stats_rmse.mean()) << std::endl;
					}
				}
				if(para.learnflag<99) write_single_ctx(single_ctxunit_files, single_output_files, time, para);
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

	gsl_odeiv2_driver * RNNdriver_test = gsl_odeiv2_driver_alloc_y_new (&RNNsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

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
			euler_step_cortex(time, hcx, para);
			// compute output
			para.Z=output_rate(para.Wout, para.nuCx);

			// compute test error
			error_signed_test = para.Z-func_target(time, para);
			RMSE_test+=error_signed_test*error_signed_test;

			// get open closed trajectory
			if( (ti < para.Tmax_FM) && bool_traj) {
				for(int i= 0; i < para.Ncx; i++) 			INPUTS_CLOSED.at(i, ti)= hcx[i];
				ti++;
			}

			// write stuff to files
			if ((t % (para.write_every* int(para.dt/para.dt_euler)))==0) {
				if(para.learnflag<99){
					error_signed_test_file << error_signed_test << std::endl;

					for(int i= 0; i < para.Ncx; i++)	nuCx_test_file << para.nuCx[i] << " ";
					nuCx_test_file << std::endl;
				}
				if(para.ROwriteflag==1) {
					Z_test_files[np] << para.Z << "\t";
				}
			}

			if(para.learnflag<99) write_single_ctx(single_ctxunit_files, single_output_files, time, para);
			tt++;
		}
		// calculate test RMSE
		RMSE_test=sqrt(RMSE_test/double(para.Tmax_test_euler))/sqrt(para.signalvar);
		alg_steps_file << "\t \t \t End of test phase, T =  "<< para.Tmax_test_euler*para.dt_euler << "\t Ttot =  "<< time << std::endl;
	}
	// RK integration (noise-free)
	else{
		for(int t= 0; t < para.Tmax_test; t++){

			status = gsl_odeiv2_driver_apply (RNNdriver_test, &time, tt*para.dt, hcx);
			update_rates_ctx(time, hcx, para);
			if (status != GSL_SUCCESS){
				printf ("Error in test phase, return value=%d\n", status);
				break;
			}

			// compute output
			para.Z=output_rate(para.Wout, para.nuCx);

			// compute test error
			error_signed_test = para.Z-func_target(time, para);
			RMSE_test+=error_signed_test*error_signed_test;

			// get closed loop trajectory
			if(t < para.Tmax_FM) {
				for(int i= 0; i < para.Ncx; i++) INPUTS_CLOSED.at(i, t)= hcx[i];
			}

			// write stuff to files
			if ((t % para.write_every)==0) {
				if(para.learnflag<99){
					error_signed_test_file << error_signed_test << std::endl;

					for(int i= 0; i < para.Ncx; i++)	nuCx_test_file << para.nuCx[i] << " ";
					nuCx_test_file << std::endl;
				}
				if(para.ROwriteflag==1) {
					Z_test_files[np] << para.Z << "\t";
				}
			}

			if(para.learnflag<99) write_single_ctx(single_ctxunit_files, single_output_files, time, para);
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
		Z_test_files[np] << std::endl;
	}

	gsl_odeiv2_driver_free (RNNdriver_test);

//	**************************************************
// 	Calculate Floquet multipliers
//	**************************************************
	//~ calc_floquet_RK(hcx_save, para, json_measures, alg_steps_file);

	calc_floquet_RK_vector(hcx_save, para, json_measures, alg_steps_file);

//	*************************
// 	Calculate norm of output weights
//	*************************
	double Wnorm_ctx=calc_wnorm(para.Wout);

//	**************************************************
// 	Calculate Fourier components
//	**************************************************
	gsl_odeiv2_driver * RNNdriver_FT = gsl_odeiv2_driver_alloc_y_new (&RNNsys, gsl_odeiv2_step_rkf45,para.dt, para.INT_delta_abs,para.INT_delta_rel);

	// set loop parameter
	para.closedloop=0;		// open loop
	para.testphaseflag=0;	// no noise
	para.RKflag=1;			// Runke-Kutta integration
	alg_steps_file << "\t \t \t Transient phase for FT calculations"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	// set initial network state
	for(int i= 0; i < para.Ncx; i++) {
		hcx[i]=hcx_save[i];
	}
	stats_hcx.reset(); 	stats_nuCx.reset();
	time=0;

	tt=para.Tmax_trans;
	status = gsl_odeiv2_driver_apply (RNNdriver_FT, &time, (para.Tmax_trans-1)*para.dt, hcx);
	if (status != GSL_SUCCESS){
		printf ("Error during Fourier transient phase, return value=%d\n", status);
		break;
	}
	alg_steps_file << "\t \t \t End of FT transient phase, T =  "<< para.Tmax_trans*para.dt << "\t Ttot =  "<< time << std::endl;

	alg_steps_file << "\t \t \t FT trajectory extraction phase"<< "\t RKflag = "<< para.RKflag <<  "\t testphaseflag = "<< para.testphaseflag << " \t closedloop = "<< para.closedloop << std::endl;
	// multiple periods for FT trajectory
	for(int t= 0; t < para.Tmax_FT; t++){
		nuCx_avrg_pop[t]=0;

		status = gsl_odeiv2_driver_apply (RNNdriver_FT, &time, tt*para.dt, hcx);
		update_rates_ctx(time, hcx, para);
		if (status != GSL_SUCCESS){
			printf ("Error during FT phase, return value=%d\n", status);
			break;
		}

		// get data for FT
		for(int i= 0; i < para.Ncx; i++) RATES_FT.at(i, t)= para.nuCx[i];

		// get open loop trajectory
		if(t < para.Tmax_FM) {
			for(int i= 0; i < para.Ncx; i++) INPUTS_OPEN.at(i, t)= hcx[i];
		}

		// set activity of output equal to target
		para.Z =func_target(time, para);

		for(int i= 0; i < para.Ncx; i++){
			stats_hcx(hcx[i]);
			stats_nuCx(para.nuCx[i]);
			nuCx_avrg_pop[t]+=para.nuCx[i];
		}

		// write to files
		if(para.learnflag<99){
			if ((t % para.write_every)==0) {
				for(int i= 0; i < para.Ncx; i++)	nuCx_FT_file << para.nuCx[i] << " ";
				nuCx_FT_file << std::endl;
			}
		}

		if(para.learnflag<99) write_single_ctx(single_ctxunit_FT_files, single_output_FT_files, time, para);
		tt++;
	}
	alg_steps_file << "\t \t \t End of FT trajectory extraction phase, T =  "<< para.Tmax_FT*para.dt  << "\t Ttot =  "<< time << std::endl;

	gsl_odeiv2_driver_free (RNNdriver_FT);

//	*************************
// 	Calculate measures of interest
//	*************************
	// calculate fluctuation of population mean
	nuCx_avrg_pop=nuCx_avrg_pop/double(para.Ncx);
	double nuCx_avrg_std=arma::stddev(nuCx_avrg_pop);
	double nuCx_avrg_cv=nuCx_avrg_std/arma::mean(nuCx_avrg_pop);

	// write stats to file
	std::ofstream stats_file(para.folder+"stats.dat");		if (!stats_file) exit(EXIT_FAILURE);

	json_stats["hcx_mean"].push_back(stats_hcx.mean());
	json_stats["hcx_std"].push_back(stats_hcx.stddev());
	json_stats["nucx_mean"].push_back(stats_nuCx.mean());
	json_stats["nucx_std"].push_back(stats_nuCx.stddev());
	json_stats["nucx_pop_std"].push_back(nuCx_avrg_std);
	json_stats["nucx_pop_cv"].push_back(nuCx_avrg_cv);

	stats_file << std::setw(4) << json_stats << '\n';
	stats_file.close();

	// calculate Fourier components & FT dimensionality
	alg_steps_file << "\t \t \t Fourier analysis" << std::endl;
	FTcomp FTcomp_ctx;
	calc_FT(RATES_FT, FTcomp_ctx, para, para.Tmax_FT,"ctx");//,para.sigmaLS*para.FACsigmaLS_ctx);		// ADJUST LAST ARGUMENT FOR RLS CASE (i.e. test-phase noise)

	// calculate dimensionality
	alg_steps_file << "\t \t \t Calculate dimensionality" << std::endl;
	double DIM_ctx=calc_dim_noEV(RATES_FT,para.Tmax_FT);

	// write to files: ctx
	alg_steps_file << "\t \t \t Write population measures to files"<< std::endl;
	json_measures["DIM_ctx"].push_back(DIM_ctx);
	json_measures["DIM_FT_ctx"].push_back(FTcomp_ctx.dim_signal);
	json_measures["COND_ctx"].push_back(FTcomp_ctx.cond_signal);
	json_measures["COND_total_ctx"].push_back(FTcomp_ctx.cond_total);

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
		std::ofstream Wcx_file(para.path_connect+"Wcx.dat");					if (!Wcx_file) exit(EXIT_FAILURE);
		std::ofstream Wout_file(para.path_connect+"Wout.dat");					if (!Wout_file) exit(EXIT_FAILURE);
		std::ofstream Wfb_file(para.path_connect+"Wfb.dat");					if (!Wfb_file) exit(EXIT_FAILURE);
		std::ofstream Wfb_interf_file(para.path_connect+"Wfb_interf.dat");		if (!Wfb_interf_file) exit(EXIT_FAILURE);

		//	write stuff to connectivity-files
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
		}
		for(int i= 0; i < para.Ncx; i++) {
			json_trial["tau_distrb"].push_back(para.tau_distrb[i]);
		}

		FTcomp_ctx.Angles_all.save(angles_ctx_file,arma::raw_ascii);
		FTcomp_ctx.S_total.save(S_total_ctx_file,arma::raw_ascii);

		// write to file
		std::ofstream single_trial_file(para.folder+"single_trial.dat");		if (!single_trial_file) exit(EXIT_FAILURE);
		single_trial_file << std::setw(4) << json_trial << '\n';
		single_trial_file.close();

		Wcx_file.close();
		Wout_file.close();
		Wfb_file.close();
		Wfb_interf_file.close();
	}

	bool write_FULL=1;
	alg_steps_file << "\t \t \t Write ctx summary data to files"<< std::endl;
	write_measures_final(json_measures, FTcomp_ctx, RMSE_train, RMSE_test, Wnorm_ctx, "ctx", para, write_FULL);

	simcount_file << 1 << std::endl;

	} // end of Nreal-loop

	Z_test_files[np].close();

	} // end of Npara-loop

//	*************************
// 	terminate program
	alg_steps_file << "Terminate"<< std::endl;
   	nuCx_learn_file.close();
   	nuCx_test_file.close();
	nuCx_FT_file.close();

	Z_file.close();
	Inputs_file.close();

	error_signed_test_file.close();
	error_signed_train_file.close();
	error_full_train_file.close();
	LR_file.close();

	single_ctxunit_files.close();
	single_output_files.close();
	single_ctxunit_FT_files.close();
	single_output_FT_files.close();

	freq_file.close();
	amp_file.close();
	phase_file.close();
	freq_intf_file.close();
	amp_intf_file.close();
	phase_intf_file.close();
	angles_ctx_file.close();
	S_total_ctx_file.close();

	alg_steps_file.close();
	simcount_file.close();

	gsl_rng_free (para.r);

  return(0);

}
