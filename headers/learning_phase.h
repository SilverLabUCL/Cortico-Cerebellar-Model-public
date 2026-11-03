void learning_phase(arma::colvec& GC, signal &ft, struct RLS &RLS_cereb, presyn J[], double PC[], double MAE_train[], double RMSE_train[], double error_signed_train[], netpara &para) {


(arma::colvec& nuCxm, arma::colvec& nuCx, arma::mat& inputsignal, double hcx[], presyn Wcx[], presyn Wfb[], presyn Wout[], double Z[], double PC[], gsl_rng * r, int &t, netpara &para, int &algflag, int &loopflag, int &noiseflag)


	int i,k;
	int Ncx=para.Ncx;
	int Ngc=para.Ngc;
	int Nmf=para.Nmf;

//  reset initial conditions
	for(i= 0; i < Ncx; i++) {
		hcx[i]=hcxini[i];
		nuCx[i]=(*para.Ptr_phi)(hcx[i]);
	}
	for(i= 0; i < Nmf; i++) MF[i]=nuCx[Ncx_to_MF_idx[i]];
	for(i= 0; i < Ngc; i++) GC[i]=hGC[i]=fsparse[i]=0;
	for(k= 0; k < Npc; k++) PC[k]=Z[k]=MAE_train[k]=RMSE_train[k]=error_signed_train[k]=error_signed_test[k]=0;
	//~ tsvd=0;
	for(t= 0; t < Tmax_learn; t++){
		// write stuff to files
		if ((t % write_every)==0 && learnflag<99) {
			for(i= 0; i < Nmf; i++) MF_file << MF[i] << " ";
			MF_file << std::endl;

			for(i= 0; i < Ncx; i++)	nuCx_learn_file << nuCx[i] << " ";
			nuCx_learn_file << std::endl;

			for(i= 0; i < Ngc; i++) GC_learn_file << GC[i] << " ";
			GC_learn_file << std::endl;

			for(k= 0; k < Npc; k++) {
				PC_file << PC[k] << " ";
				Z_file << Z[k] << " ";
			}
			PC_file << std::endl;
			Z_file  << std::endl;
		}
		euler_step_cortex(nuCxm, nuCx, ft.target, hcx, Wcx, Wfb, Wout, Z, PC, r, t, para, algflag, loopflag, noiseflag);

		// cerebellar populations
		if(learnflag>=99 && algflag==1){
		}
		else{
			euler_step_cerebellum(nuCxm, GC, hGC, MF, Ncx_to_MF_idx, theta, gain, Wmf, J, PC, r, para, noiseflag);
		}
		if ((t % learn_every)==0) {
			// FORCE/DELTA learning with cortical population
			if(algflag==1){
				if(LRflag==0)		RLS_step_cortex(nuCx, ft, RLS_ctx, Wout, Z, MAE_train, RMSE_train, error_signed_train, t, para);
				else if(LRflag==1)	DELTA_step(nuCx, ft, RLS_ctx, Wout, Z, MAE_train, RMSE_train, error_signed_train, t, para);

			}
			// FORCE/DELTA learning with cerebellar cortex
			if(algflag==2){
				if(LRflag==0)		RLS_step_cerebellum(GC, ft, RLS_cereb, J, PC, MAE_train, RMSE_train, error_signed_train, t, para);
				else if(LRflag==1)	DELTA_step(GC, ft, RLS_cereb, J, PC, MAE_train, RMSE_train, error_signed_train, t, para);
			}
		}
		if (learnflag<99) {
			if ((t % write_every)==0) {
				if(algflag==1)		{
					calc_wnorm_change(RLS_ctx, Ncx, para);
					LR_file << RLS_ctx.LR.at(0,0) << "\t" << RLS_ctx.LR.at(0,1) << "\t" << RLS_ctx.LR.at(0,2) << std::endl;
				}
				else if(algflag==2)	{
					calc_wnorm_change(RLS_cereb, Ngc, para);
					LR_file << RLS_cereb.LR.at(0,0) << "\t" << RLS_cereb.LR.at(0,1) << "\t" << RLS_cereb.LR.at(0,2) << std::endl;
				}

				for(k= 0; k < Npc; k++) error_signed_train_file << error_signed_train[k] << "\t";
				error_signed_train_file << std::endl;

				for(i= 0; i < 30; i++) {
					theta_temp_file << theta[i] << "\t";
					J_temp_file << J[0].weight[i] << "\t";
					Wmf_temp_file << Wmf[i].weight[0] << "\t";
					gain_temp_file << gain[i]<< "\t";
				}
				theta_temp_file << std::endl;
				Wmf_temp_file << std::endl;
				J_temp_file << std::endl;
				gain_temp_file << std::endl;
			}
		}
	}

	return;
}
