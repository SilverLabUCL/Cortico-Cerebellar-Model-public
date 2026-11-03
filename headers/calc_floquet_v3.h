// added Mf to GC expansion
// all weights in para struc

void calc_floquet(arma::mat& rates_ctx, arma::mat& inputs_ctx, arma::mat& inputs_gc, netpara &para, nlohmann::json &json_measures, std::ofstream &alg_steps_file) {

if(para.floquetflag==1){
	if(para.algswitch==1)	alg_steps_file << "\t \t \t Calculate Floquet multipliers for ctx"<< std::endl;
	else					alg_steps_file << "\t \t \t Calculate Floquet multipliers for ctx + cereb"<< std::endl;

	int Tsteps = round(1.0/(para.lambda*para.dt));		// one period
	//~ std::cout << "Tsteps = "<< Tsteps << std::endl;

	arma::vec dt_tau=arma::ones<arma::vec>(para.Ncx);
	if(para.taudistrbflag==1)	dt_tau=para.dttau_distrb;
	else						dt_tau=dt_tau*para.dt/para.tau_avrg;

	arma::vec dt_tauGC=arma::ones<arma::vec>(para.Ngc);
	if(para.taugcmodeflag==2)	dt_tauGC=para.dttauGC_distb;
	else						dt_tauGC=dt_tauGC*para.dt/para.tauGC_avrg;

	//~ std::cout << para.tau_avrg <<"\t "<< para.dt/para.tau_avrg << std::endl;
	//~ std::cout << para.tauGC_avrg <<"\t "<< para.dt/para.tauGC_avrg << std::endl;

	// recurrent weights
	arma::mat Wcortex= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
	for(int i= 0; i < para.Ncx; i++) {
		for (int j=0; j<para.Wcx[i].Nsyn; j++) {
			Wcortex.at(i,para.Wcx[i].idx[j])=para.Wcx[i].weight[j];
		}
	}
	arma::mat Sm;
	arma::mat Iall;
	arma::mat I = arma::eye( para.Ncx, para.Ncx );

	arma::mat dPHI(arma::size(inputs_ctx));
	arma::mat dPSI(arma::size(inputs_gc));

	//~ std:: cout << arma::size(inputs_ctx) << std::endl;

	//~ std::cout << inputs_ctx(0,0) <<"\t "<< inputs_ctx(0,Tsteps) <<"\t "<< rates_ctx(0,0) <<"\t "<< rates_ctx(0,Tsteps) << std::endl;
	//~ std:: cout << arma::norm(inputs_ctx.col(0) - inputs_ctx.col(Tsteps)) << std::endl;

	for(int i= 0; i < 1; i++) {
		if(para.algswitch==1)	json_measures["Deviation_LC_ctx"].push_back(arma::norm(inputs_ctx.col(0) - inputs_ctx.col(Tsteps-1)));
		else					json_measures["Deviation_LC_cereb"].push_back(arma::norm(inputs_gc.col(0) - inputs_gc.col(Tsteps-1)));
	}

	// tanh derivative
	if(para.actfuncflag==0)		dPHI=1.0/(arma::cosh(inputs_ctx) % arma::cosh(inputs_ctx));
	// cortex only
	//~ if(para.algswitch==1){
	if(para.cortexonlyflag==1 && para.algswitch==1){		// only compute FM when in CORTEXonly mode
		Sm = arma::eye( para.Ncx, para.Ncx );
		arma::mat Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);
		// linear RNN with interference signal
		// if(para.actfuncflag==1 && para.intrfflag==1){
		if(para.actfuncflag==1){
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					Weff.at(i,j)=-I.at(i,j)+Wcortex.at(i,j) + para.Wfb[i]*para.Wout[j];
				}
			}
			Weff= Weff*arma::diagmat(dt_tau);
			// std::cout << "muh" << std::endl;
			for(int tt= 0; tt < Tsteps; tt++) {
				//~ Sm=Sm +Weff *Sm;
				Sm=(I +Weff) *Sm;
			}
		}
		// tanh RNN
		else if(para.actfuncflag==0){
			arma::mat WrankN( para.Ncx, para.Ncx );
			for(int i= 0; i < para.Ncx; i++) {
				for (int j=0; j<para.Ncx; j++) {
					WrankN.at(i,j)=Wcortex.at(i,j) + para.Wfb[i]*para.Wout[j];
				}
			}
			for(int tt= 0; tt < Tsteps; tt++) {
				Weff=-I+WrankN*arma::diagmat(dPHI.col(tt));
				//~ Sm=Sm +Weff*arma::diagmat(dt_tau) *Sm;
				Sm=(I +Weff*arma::diagmat(dt_tau)) *Sm;
			}
		}
	}
	else {
		Sm = arma::eye( para.Ncx, para.Ncx );
	}
	// cortex + cerebellum
	//~ else{
	if(para.algswitch==2){
		// MF-GC weights
		arma::mat Wmfgc= arma::zeros<arma::mat>(para.Ngc,para.Nmf);
		for(int i= 0; i < para.Ngc; i++) {
			for (int j=0; j<para.Wmf[i].Nsyn; j++) {
				Wmfgc.at(i,para.Wmf[i].idx[j])=1.0/sqrt(para.Nsyn);
			}
		}
		// CTX-MF weights
		arma::mat Wctxmf= arma::zeros<arma::mat>(para.Nmf,para.Ncx);
		for(int i= 0; i < para.Nmf; i++) {
			Wctxmf.at(i,para.Ncx_to_MF_idx[i])=1.0;
		}
		// effective CTX-GC weights
		arma::mat Wctxgc=Wmfgc*Wctxmf;

		// rank 1 matrix
		arma::mat Wrank1( para.Ncx, para.Ngc );
		for(int i= 0; i < para.Ncx; i++) {
			for (int j=0; j<para.Ngc; j++) {
				Wrank1.at(i,j)= para.Wfb[i]*para.Wout_PC[j];
			}
		}

		// ReLu derivative (heaviside)
		if(para.actfuncflag_GC==0){
			for(int i= 0; i < para.Ngc; i++) {
				for(int tt= 0; tt < Tsteps; tt++) {
					dPSI.at(i,tt)=para.gain_vector[i]*heavi(inputs_gc.at(i,tt)-para.threshold_vector[i]);
				}
			}
		}
		// softplus derivative
		else if(para.actfuncflag_GC==1){
			for(int i= 0; i < para.Ngc; i++) {
				for(int tt= 0; tt < Tsteps; tt++) {
					dPSI.at(i,tt)=para.gain_vector[i]*inv_logit(inputs_gc.at(i,tt)-para.threshold_vector[i], para.sp_alpha);
				}
			}
		}
		// *************************************
		// instantaneous GCs
		// *************************************
		if(para.taugcmodeflag==0) {
			Sm = arma::eye( para.Ncx, para.Ncx );
			arma::mat Weff= arma::zeros<arma::mat>(para.Ncx,para.Ncx);

			// rank N matrix
			arma::mat WrankN=(-I+Wcortex)*arma::diagmat(dt_tau);

			// without Golgi cell feedback
			if(para.gocmodeflag==0){
				if(para.actfuncflag==1 && para.actfuncflag_GC==2){
					// linear RNN & cereb linear
					Wrank1=Wrank1*arma::diagmat(para.gain_vector)*Wctxgc*arma::diagmat(dt_tau);
					Weff=WrankN+Wrank1;
					for(int tt= 0; tt < Tsteps; tt++) {
						//~ Sm=Sm +Weff *Sm;
						Sm=(I +Weff) *Sm;
					}
				}
				else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					Wctxgc=Wctxgc*arma::diagmat(dt_tau);
					for(int tt= 0; tt < Tsteps; tt++) {
						Weff=WrankN+Wrank1*arma::diagmat(dPSI.col(tt))*Wctxgc;
						//~ Sm=Sm +Weff*Sm;
						Sm=(I +Weff) *Sm;
					}
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						Weff=-I+(Wcortex+Wrank1*arma::diagmat(dPSI.col(tt))*Wctxgc)*arma::diagmat(dPHI.col(tt));
						//~ Sm=Sm +Weff*arma::diagmat(dt_tau) *Sm;
						Sm=(I +Weff*arma::diagmat(dt_tau)) *Sm;
					}
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
			// with Golgi cell feedback
			else if (para.gocmodeflag==1){
				arma::mat Igc= arma::eye( para.Ngc, para.Ngc );
				arma::mat Wtmp( para.Ngc, para.Ngc );
				arma::mat Wgoc_distrb_vec=para.Wgoc_distrb/double(para.Ngc);
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					Wctxgc=Wctxgc*arma::diagmat(dt_tau);
					for(int tt= 0; tt < Tsteps; tt++) {
						Wtmp=arma::inv(Igc+Wgoc_distrb_vec*arma::trans(dPSI.col(tt)));
						Weff=WrankN+Wrank1*arma::diagmat(dPSI.col(tt))*Wtmp*Wctxgc;
						//~ Sm=Sm +Weff*Sm;
						Sm=(I +Weff) *Sm;
					}
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						//~ Wtmp=arma::inv(Igc+Wgoc_distrb_mat*arma::repmat(dPSI.col(tt),1,para.Ngc));
						Wtmp=arma::inv(Igc+Wgoc_distrb_vec*arma::trans(dPSI.col(tt)));
						Weff=-I+(Wcortex+Wrank1*arma::diagmat(dPSI.col(tt))*Wtmp*Wctxgc)*arma::diagmat(dPHI.col(tt));
						//~ Sm=Sm +Weff*arma::diagmat(dt_tau)*Sm;
						Sm=(I +Weff*arma::diagmat(dt_tau)) *Sm;
					}
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
			else if (para.gocmodeflag==2)	std::cerr << "Heterogeneous Golgi cell feedback not implemented (yet)." << std::endl;
			// with feedforward Golgi cell inhibition
			else if (para.gocmodeflag==3)	{
				arma::vec SD_rate_mfs=arma::stddev(rates_ctx.t(),0,1);
				arma::vec avrg_rate_mfs=arma::mean(rates_ctx.t(),1);
				arma::mat Wgoc_distrb_vec=para.Wgoc_distrb/double(para.Ngc);
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					WrankN=(-I+Wcortex);
					for(int tt= 0; tt < Tsteps; tt++) {
						//~ Wtmp = Wgoc_distrb_vec*arma::trans(rates_ctx.col(tt)-avrg_rate_mfs[tt])/(SD_rate_mfs[tt]);
						//~ Weff=WrankN +Wrank1*arma::diagmat(dPSI.col(tt))*(Wctxgc -Wtmp);
						Weff=WrankN +Wrank1*arma::diagmat(dPSI.col(tt))*(Wctxgc -Wgoc_distrb_vec*arma::trans(rates_ctx.col(tt)-avrg_rate_mfs[tt])/SD_rate_mfs[tt]);
						//~ Sm=Sm +Weff*arma::diagmat(dt_tau)*Sm;
						Sm=(I +Weff*arma::diagmat(dt_tau)) *Sm;
					}
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						//~ Wtmp = Wgoc_distrb_vec*arma::trans(dPHI.col(tt)%(rates_ctx.col(tt)-avrg_rate_mfs[tt]))/(SD_rate_mfs[tt]);
						//~ Weff=-I+Wcortex*arma::diagmat(dPHI.col(tt))+Wrank1*arma::diagmat(dPSI.col(tt))*(Wctxgc*arma::diagmat(dPHI.col(tt)) -Wtmp);
						Weff=-I+Wcortex*arma::diagmat(dPHI.col(tt))+Wrank1*arma::diagmat(dPSI.col(tt))*(Wctxgc*arma::diagmat(dPHI.col(tt)) -Wgoc_distrb_vec*arma::trans( dPHI.col(tt)%(rates_ctx.col(tt)-avrg_rate_mfs[tt]))/(SD_rate_mfs[tt]));
						//~ Sm=Sm +Weff*arma::diagmat(dt_tau)*Sm;
						Sm=(I +Weff*arma::diagmat(dt_tau)) *Sm;
					}
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
		}
		// *************************************
		// non-instantaneous GCs
		// *************************************
		else{
			Sm = arma::eye( para.Ncx + para.Ngc, para.Ncx + para.Ngc);
			Iall = arma::eye( para.Ncx + para.Ngc, para.Ncx + para.Ngc);
			arma::mat Weff( para.Ncx + para.Ngc, para.Ncx + para.Ngc);
			arma::vec dttauALL = arma::join_vert(dt_tau,dt_tauGC);
			// rank N matrix
			arma::mat WrankN=-I + Wcortex;
			arma::mat Igc= arma::eye( para.Ngc, para.Ngc );
			arma::mat AB;
			arma::mat CD;

			// without Golgi cell feedback
			if(para.gocmodeflag==0){
				if(para.actfuncflag==1 && para.actfuncflag_GC==2){
					// linear RNN & cereb linear
					AB = arma::join_horiz(WrankN,Wrank1);
					CD = arma::join_horiz(Wctxgc,-Igc);
					Weff = arma::join_vert(AB,CD);
					for(int tt= 0; tt < Tsteps; tt++) {
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL))*Sm;
					}
				}
				else if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					CD = arma::join_horiz(Wctxgc,-Igc);
					for(int tt= 0; tt < Tsteps; tt++) {
						AB = arma::join_horiz( WrankN, Wrank1*arma::diagmat(dPSI.col(tt)));
						Weff=arma::join_vert( AB, CD);
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL))*Sm;
					}
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						AB = arma::join_horiz( -I + Wcortex*arma::diagmat(dPHI.col(tt)), Wrank1*arma::diagmat(dPSI.col(tt)));
						CD = arma::join_horiz(Wctxgc*arma::diagmat(dPHI.col(tt)),-Igc);
						Weff=arma::join_vert( AB, CD);
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL)) *Sm;
					}
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}
			// with Golgi cell feedback
			else if (para.gocmodeflag==1){
				arma::mat Wgoc_distrb_vec=para.Wgoc_distrb/double(para.Ngc);
				arma::mat Wtmp= arma::zeros<arma::mat>(para.Ngc,para.Ngc);
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						Wtmp = Igc+Wgoc_distrb_vec*arma::trans(dPSI.col(tt));
						AB = arma::join_horiz( WrankN, Wrank1*arma::diagmat(dPSI.col(tt)));
						CD = arma::join_horiz(Wctxgc,-Wtmp);
						Weff=arma::join_vert( AB, CD);
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL)) *Sm;
					}
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						Wtmp = Igc+Wgoc_distrb_vec*arma::trans(dPSI.col(tt));
						AB = arma::join_horiz( -I+Wcortex*arma::diagmat(dPHI.col(tt)), Wrank1*arma::diagmat(dPSI.col(tt)));
						CD = arma::join_horiz(Wctxgc*arma::diagmat(dPHI.col(tt)),-Wtmp);
						Weff=arma::join_vert( AB, CD);
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL)) *Sm;
					}
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;

			}
			else if (para.gocmodeflag==2)	std::cerr << "Heterogeneous Golgi cell feedback not implemented (yet)." << std::endl;
			else if (para.gocmodeflag==3){
				arma::vec SD_rate_mfs=arma::stddev(rates_ctx.t(),0,1);
				arma::vec avrg_rate_mfs=arma::mean(rates_ctx.t(),1);
				arma::mat Wgoc_distrb_vec=para.Wgoc_distrb/double(para.Ngc);
				arma::mat Wtmp= arma::zeros<arma::mat>(para.Ngc,para.Ncx);
				if(para.actfuncflag==1 && para.actfuncflag_GC<2){
					// linear RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						Wtmp = Wgoc_distrb_vec*arma::trans(rates_ctx.col(tt)-avrg_rate_mfs[tt])/(SD_rate_mfs[tt]);
						AB = arma::join_horiz( WrankN, Wrank1*arma::diagmat(dPSI.col(tt)));
						CD = arma::join_horiz(Wctxgc-Wtmp,-Igc);
						Weff=arma::join_vert( AB, CD);
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL)) *Sm;
					}
				}
				else if(para.actfuncflag==0 && para.actfuncflag_GC<2){
					// tanh RNN & cereb ReLu or softplus
					for(int tt= 0; tt < Tsteps; tt++) {
						Wtmp = Wgoc_distrb_vec*arma::trans( dPHI.col(tt)%(rates_ctx.col(tt)-avrg_rate_mfs[tt]))/(SD_rate_mfs[tt]);
						AB = arma::join_horiz( -I+Wcortex*arma::diagmat(dPHI.col(tt)), Wrank1*arma::diagmat(dPSI.col(tt)));
						CD = arma::join_horiz(Wctxgc*arma::diagmat(dPHI.col(tt))-Wtmp,-Igc);
						Weff=arma::join_vert( AB, CD);
						//~ Sm=Sm +Weff*arma::diagmat(dttauALL)*Sm;
						Sm=(Iall +Weff*arma::diagmat(dttauALL)) *Sm;
					}
				}
				else 	std::cerr << "Combination of input-output relations not implemented (yet)." << std::endl;
			}

		}
	}

	//~ for(int i= 0; i < 5; i++) {
		//~ for (int j=0; j<5; j++) {
			 //~ std::cout << Sm.at(i,j) << "\t"<<std::endl;
			//~ // std::cout << dPHI.at(i,j) << "\t"<<std::endl;
			//~ // std::cout << dPSI.at(i,j) << "\t"<<std::endl;
		//~ }
		//~ std::cout << std::endl;
	//~ }
	//~ std::cout << Sm.has_nan() << std::endl;

	// calculate EVs of monodromy matrix Sm
	int NEV=para.Ncx;
	if(para.taugcmodeflag>0) NEV=para.Ncx+para.Ngc;
	arma::cx_vec Sm_eigval(NEV);
	if(Sm.has_nan()==0){
		Sm_eigval = arma::eig_gen( Sm );
		Sm_eigval=Sm_eigval(arma::sort_index( arma::abs(Sm_eigval), "descend" ));
	}
	else 	Sm_eigval(0)=arma::datum::nan;
	//~ std::cout << arma::size(Sm_eigval) << std::endl;

	//~ if(para.algswitch==1)	std::cout << "RNN Floquet multipliers:" << std::endl;
	//~ else 					std::cout << "Cereb Floquet multipliers:" << std::endl;
	//~ for(int i= 0; i < 10; i++) std::cout << std::real(Sm_eigval[i]) <<  "\t" << std::imag(Sm_eigval[i]) <<  "\t" << std::abs(Sm_eigval[i]) << std::endl;
	//~ std::cout << std::endl;
	for(int i= 0; i < 1; i++) {
		if(para.algswitch==1)	json_measures["Max_Floquet_multiplier_ctx"].push_back(std::abs(Sm_eigval[i]));
		else					json_measures["Max_Floquet_multiplier_cereb"].push_back(std::abs(Sm_eigval[i]));
	}

}
else{
	if(para.algswitch==1)	json_measures["Max_Floquet_multiplier_ctx"].push_back(NAN);
	else					json_measures["Max_Floquet_multiplier_cereb"].push_back(NAN);
}

	return;
}
