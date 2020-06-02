//////////////////////////////////////////////////////////
///// Structure for saveing post samplings /////
//////////////////////////////////////////////////////////

struct Sampling{
    
    double ** basis2;
    double ** zb;
    double ** zb2;
    double ** xb;
    double ** xb_test;
    double ** rxb;
    
    // Intercept for outcomes
    double ** Zinter;
    double ** mean_Zinter;

    // Theta
    double ** theta;
    double ** mean_theta;
    
    // Loading
    double ** load;
    double ** loadstar;
    double ** mean_load;
    double ** mean_loadstar;
    
    // Latent
    double ** latent;
    double ** latentstar;
    double ** latent2;
    double ** mean_latent;
    double ** mean_latentstar;
    
    // Phi
    double * Phi2Inv;
    double * mean_Phi2Inv;
    double Phi2Inv_a;
    double Phi2Inv_b;
    
    // Beta
    double ** beta;
    double ** betastar;
    double ** mean_beta;
    double ** mean_betastar;
    
    // Alpha
    double ** alpha;
    double ** alphastar;
    double ** mean_alpha;
    double ** mean_alphastar;
    
    // XBA, XBAR
    double ** xba; // X*Baisi*Aalpha
    double ** xbar; // X*Baisi*Aalpha*Gamma
    
    // Mu
    double ** mustar;
    double ** mean_mu;
    
    // Summarized image
    double ** sumX;
    
    // Predictor parameters
    double * gamma;
    double * mean_gamma;
    double * mean_gamma_update;
    double * w;
    double * mean_w;
    double w_a;
    double w_b;
    
    // Errors
    double * err2inv_e;
    double * err2inv_e_mean;
    double * err2inv_zeta;
    double * err2inv_zeta_mean;
    double * err2inv_eps;
    double * err2inv_eps_mean;
    double * err2inv_u;
    double * err2inv_u_mean;
    
    // prior of errors
    double err2inv_e_a;
    double err2inv_e_b;
    double err2inv_zeta_a;
    double err2inv_zeta_b;
    double err2inv_eps_a;
    double err2inv_eps_b;
    double err2inv_alpha;
    double err2inv_load;
    double err2inv_mu;
    double err2inv_u_a;
    double err2inv_u_b;
    
    // Testing Set
    double ** Zinter_test;
    double ** latent_test;
    double ** theta_test;
    double ** outcome_test;
    double ** Zinter_test_mean;
    double ** latent_test_mean;
    double ** theta_test_mean;
    double ** outcome_test_mean;
    double ** outcome_test_mean2;
    
    // Training Set
    double ** latent_train;
    double ** theta_train;
    double ** outcome_train;
    double ** latent_train_mean;
    double ** theta_train_mean;
    double ** outcome_train_mean;
    double ** outcome_train_mean2;
    double ** latent_train_err;
    double ** theta_train_err;
};


struct Sampling setupSamp(const int M, const int nobs, const int nts, const int L, const int P, const int K, struct BasisFunc BF, struct Inputdata data){
    
    struct Sampling PostSamp;
    int l, Ml, parcel_len;
    
    ///// Allocate Memory /////
    // data
    PostSamp.basis2 = (double **)malloc(L * sizeof(double *));
    PostSamp.zb = (double **)malloc(L * sizeof(double *));
    PostSamp.zb2 = (double **)malloc(L * sizeof(double *));
    PostSamp.xb = (double **)malloc(L * sizeof(double *));
    PostSamp.rxb = (double **)malloc(L * sizeof(double *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        PostSamp.basis2[l] = (double *)malloc(Ml * Ml * sizeof(double));
        PostSamp.zb[l] = (double *)malloc(nobs * Ml * sizeof(double));
        PostSamp.zb2[l] = (double *)malloc(nobs * Ml * sizeof(double));
        PostSamp.xb[l] = (double *)malloc(nobs * P * Ml * sizeof(double));
        PostSamp.rxb[l] = (double *)malloc(nobs * Ml * sizeof(double));
    }
    if(nts > 0){
        PostSamp.xb_test = (double **)malloc(L * sizeof(double *));
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            PostSamp.xb_test[l] = (double *)malloc(nts*P*Ml* sizeof(double));
        }
    }
    
    
    
    
    // Parameter
    PostSamp.Zinter = (double **)malloc(L * sizeof(double *));
    PostSamp.theta = (double **)malloc(L * sizeof(double *));
    PostSamp.load = (double **)malloc(L * sizeof(double *));
    PostSamp.loadstar = (double **)malloc(L * sizeof(double *));
    PostSamp.latent = (double **)malloc(L * sizeof(double *));
    PostSamp.latentstar = (double **)malloc(L * sizeof(double *));
    PostSamp.latent2 = (double **)malloc(L * sizeof(double *));
    PostSamp.beta = (double **)malloc(L * sizeof(double *));
    PostSamp.betastar = (double **)malloc(L * sizeof(double *));
    PostSamp.alpha = (double **)malloc(L * sizeof(double *));
    PostSamp.alphastar = (double **)malloc(L * sizeof(double *));
    PostSamp.xba = (double **)malloc(L * sizeof(double *));
    PostSamp.xbar = (double **)malloc(L * sizeof(double *));
    PostSamp.mustar = (double **)malloc(L * sizeof(double *));
    PostSamp.sumX = (double **)malloc(L * sizeof(double *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        PostSamp.Zinter[l] = (double *)malloc(parcel_len * sizeof(double));
        PostSamp.theta[l] = (double *)malloc(Ml * nobs * sizeof(double));
        PostSamp.load[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.loadstar[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.latent[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.latentstar[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.latent2[l] = (double *)malloc(K * K * sizeof(double));
        PostSamp.beta[l] = (double *)malloc(parcel_len * K * sizeof(double));
        PostSamp.betastar[l] = (double *)malloc(parcel_len * K * sizeof(double));
        PostSamp.alpha[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.alphastar[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.xba[l] = (double *)malloc(nobs*K*P * sizeof(double));
        PostSamp.xbar[l] = (double *)malloc(nobs*K * sizeof(double));
        PostSamp.mustar[l] = (double *)malloc(nobs*K * sizeof(double));
        PostSamp.sumX[l] = (double *)malloc(nobs*parcel_len * sizeof(double));
    }
    PostSamp.Phi2Inv = (double *)calloc(K*L, sizeof(double));
    PostSamp.gamma = (double *)calloc(P*L, sizeof(double));
    PostSamp.w = (double *)calloc(L, sizeof(double));
    PostSamp.err2inv_e = (double *)calloc(1, sizeof(double));
    PostSamp.err2inv_zeta = (double *)calloc(1, sizeof(double));
    PostSamp.err2inv_eps = (double *)calloc(1, sizeof(double));
    PostSamp.err2inv_u = (double *)calloc(L, sizeof(double));
    
    
    // Posterior Mean
    PostSamp.mean_Zinter = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_theta = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_load = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_loadstar = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_latent = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_latentstar = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_beta = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_betastar = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_alpha = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_alphastar = (double **)malloc(L * sizeof(double *));
    PostSamp.mean_mu = (double **)malloc(L * sizeof(double *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        PostSamp.mean_Zinter[l] = (double *)malloc(parcel_len * sizeof(double));
        PostSamp.mean_theta[l] = (double *)malloc(Ml * nobs * sizeof(double));
        PostSamp.mean_load[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.mean_loadstar[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.mean_latent[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.mean_latentstar[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.mean_beta[l] = (double *)malloc(parcel_len * K * sizeof(double));
        PostSamp.mean_betastar[l] = (double *)malloc(parcel_len * K * sizeof(double));
        PostSamp.mean_alpha[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.mean_alphastar[l] = (double *)malloc(Ml * K * sizeof(double));
        PostSamp.mean_mu[l] = (double *)malloc(nobs * K * sizeof(double));
    }
    PostSamp.mean_Phi2Inv = (double *)calloc(K*L, sizeof(double));
    PostSamp.mean_gamma = (double *)calloc(P*L, sizeof(double));
    PostSamp.mean_gamma_update = (double *)calloc(P*L, sizeof(double));
    PostSamp.mean_w = (double *)calloc(L, sizeof(double));
    PostSamp.err2inv_e_mean = (double *)calloc(1, sizeof(double));
    PostSamp.err2inv_zeta_mean = (double *)calloc(1, sizeof(double));
    PostSamp.err2inv_eps_mean = (double *)calloc(1, sizeof(double));
    PostSamp.err2inv_u_mean = (double *)calloc(L, sizeof(double));
    
    ///// Priors
    PostSamp.err2inv_alpha = 1.0f;
    PostSamp.err2inv_load = 1.0f;
    PostSamp.err2inv_mu = 1.0f;
    PostSamp.Phi2Inv_a = 1.0f;
    PostSamp.Phi2Inv_b = 1.0f;
    PostSamp.err2inv_e_a = 1.0f;
    PostSamp.err2inv_e_b = 1.0f;
    PostSamp.err2inv_zeta_a = 1.0f;
    PostSamp.err2inv_zeta_b = 1.0f;
    PostSamp.err2inv_eps_a = 1.0f;
    PostSamp.err2inv_eps_b = 1.0f;
    PostSamp.err2inv_u_a = 1.0f;
    PostSamp.err2inv_u_b = 1.0f;
    PostSamp.w_a = 1.0f;
    PostSamp.w_b = 1.0f;
    
    
    ///// Testing
    if(nts > 0){
        PostSamp.Zinter_test = (double **)malloc(L * sizeof(double *));
        PostSamp.latent_test = (double **)malloc(L * sizeof(double *));
        PostSamp.theta_test = (double **)malloc(L * sizeof(double *));
        PostSamp.outcome_test = (double **)malloc(L * sizeof(double *));
        PostSamp.Zinter_test_mean = (double **)malloc(L * sizeof(double *));
        PostSamp.latent_test_mean = (double **)malloc(L * sizeof(double *));
        PostSamp.theta_test_mean = (double **)malloc(L * sizeof(double *));
        PostSamp.outcome_test_mean = (double **)malloc(L * sizeof(double *));
        PostSamp.outcome_test_mean2 = (double **)malloc(L * sizeof(double *));
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            PostSamp.Zinter_test[l] = (double *)malloc(parcel_len * sizeof(double));
            PostSamp.latent_test[l] = (double *)malloc(nts * K * sizeof(double));
            PostSamp.theta_test[l] = (double *)malloc(Ml * nts * sizeof(double));
            PostSamp.outcome_test[l] = (double *)malloc(parcel_len * nts * sizeof(double));
            PostSamp.Zinter_test_mean[l] = (double *)malloc(parcel_len * sizeof(double));
            PostSamp.latent_test_mean[l] = (double *)malloc(nts * K * sizeof(double));
            PostSamp.theta_test_mean[l] = (double *)malloc(Ml * nts * sizeof(double));
            PostSamp.outcome_test_mean[l] = (double *)malloc(parcel_len * nts * sizeof(double));
            PostSamp.outcome_test_mean2[l] = (double *)malloc(parcel_len * nts * sizeof(double));
        }
    }
    
    
    ///// Training
    PostSamp.latent_train = (double **)malloc(L * sizeof(double *));
    PostSamp.theta_train = (double **)malloc(L * sizeof(double *));
    PostSamp.outcome_train = (double **)malloc(L * sizeof(double *));
    PostSamp.latent_train_mean = (double **)malloc(L * sizeof(double *));
    PostSamp.theta_train_mean = (double **)malloc(L * sizeof(double *));
    PostSamp.outcome_train_mean = (double **)malloc(L * sizeof(double *));
    PostSamp.outcome_train_mean2 = (double **)malloc(L * sizeof(double *));
    PostSamp.latent_train_err = (double **)malloc(L * sizeof(double *));
    PostSamp.theta_train_err = (double **)malloc(L * sizeof(double *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        PostSamp.latent_train[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.theta_train[l] = (double *)malloc(Ml * nobs * sizeof(double));
        PostSamp.outcome_train[l] = (double *)malloc(parcel_len * nobs * sizeof(double));
        PostSamp.latent_train_mean[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.theta_train_mean[l] = (double *)malloc(Ml * nobs * sizeof(double));
        PostSamp.outcome_train_mean[l] = (double *)malloc(parcel_len * nobs * sizeof(double));
        PostSamp.outcome_train_mean2[l] = (double *)malloc(parcel_len * nobs * sizeof(double));
        PostSamp.latent_train_err[l] = (double *)malloc(nobs * K * sizeof(double));
        PostSamp.theta_train_err[l] = (double *)malloc(Ml * nobs * sizeof(double));
    }
    
    // Compute basis^T*basis
    int i, j, k, p;
    double tempsum;
    
    for(l=0; l<L; l++){
        //printf("l=%d\n", l);
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        
        for(i=0; i<Ml; i++){
            for(j=0; j<Ml; j++){
                tempsum = 0.0;
                for(k=0; k<parcel_len; k++){
                    tempsum += (double)BF.basis[l][k*Ml+i] * (double)BF.basis[l][k*Ml+j];
                }
                PostSamp.basis2[l][i*Ml+j] = tempsum;
            }
        }
    }
    
    
    // Compute z_parcel*basis
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<nobs; i++){
            for(j=0; j<Ml; j++){
                tempsum=0.0;
                for(k=0; k<parcel_len; k++){
                    tempsum += (double)data.Zl[l][i*parcel_len+k] * (double)BF.basis[l][k*Ml+j];
                }
                PostSamp.zb[l][i*Ml+j] = tempsum;
            }
        }
    }
    
    // Compute basis * X_{ip}
    for(i=0; i<nobs; i++){

        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            
            for(p=0; p<P; p++){
                for(k=0; k<Ml; k++){
                    tempsum = 0.0;
                    for(j=0; j<parcel_len; j++){
                        tempsum += (double)data.Xl[l][i*parcel_len*P+p*parcel_len+j] * (double)BF.basis[l][j*Ml+k];
                    }
                    PostSamp.xb[l][i*P*Ml + p*Ml+ k] = tempsum;
                    
                }
            }
        }
    }
    
    // Compute basis * X_test
    if(nts > 0){
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                parcel_len = data.parcel_len[l];
                
                for(p=0; p<P; p++){
                    for(k=0; k<Ml; k++){
                        
                        tempsum = 0.0;
                        for(j=0; j<parcel_len; j++){
                            tempsum += (double)data.Xl_test[l][i*parcel_len*P+p*parcel_len+j] * (double)BF.basis[l][j*Ml+k];
                        }
                        PostSamp.xb_test[l][i*P*Ml + p*Ml+ k] = tempsum;
                    }
                }
            }
        }
    }
    

    return PostSamp;
}
