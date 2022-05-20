library(rstan)
library(bayesplot)
library(ggplot2)
library(doParallel)
library(devtools)
Sys.setenv(PATH = paste("C:/rtools40/mingw64/bin/", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
registerDoParallel(4)
getDoParWorkers()

getpara<-function(filename=NULL,names.mk=NULL){
          data.mk<-read.csv(filename)
          names.study<-unique(fixlen.str(data.mk$Study))
          hash.study<-1:length(names.study)
          names(hash.study)<-names.study
          hash.mk<-1:length(names.mk)
          names(hash.mk)<-names.mk
          para<-data.frame(study.names=fixlen.str(data.mk$Study),
          Study=array(hash.study[as.character(fixlen.str(data.mk$Study))]),
                      TP=round(data.mk$Case*data.mk$Sensitivity),
                      TN=round(data.mk$Control*data.mk$Specificity),
                      Dis=data.mk$Case,
                      NDis=data.mk$Control,
                      Test=array(hash.mk[as.character(data.mk$Marker)]))
    return(list(para=para,data=data.mk))
  }


get.ub.95<-function(x){
            y=sort(x)
            len=length(x)
            return(y[round(len*0.975)])
  }

get.lb.95<-function(x){
            y=sort(x)
            len=length(x)
            return(y[round(len*0.025)])
  }

get.nma<-function(data.ls=NULL,mk.names=NULL,n.study=1,conf=95){
	        data.nma<-data.ls$data.nma
	        func.ub<-NULL
        	func.lb<-NULL
	        if(conf==95){
		                    func.ub<-get.ub.95
		                    func.lb<-get.lb.95
	        } else if(conf==90){
		                           func.ub<-get.ub.90
		                           func.lb<-get.lb.90
	        }
	        data.nma.muse<-data.nma[paste('MU[',1,',',1:length(data.ls$mk),']',sep='')]
	        data.nma.musp<-data.nma[paste('MU[',2,',',1:length(data.ls$mk),']',sep='')]
	        data.nma.rrse<-data.nma[paste('RR[',1,',',1:length(data.ls$mk),']',sep='')]
	        data.nma.rrsp<-data.nma[paste('RR[',2,',',1:length(data.ls$mk),']',sep='')]
	        data.nma.orse<-data.nma[paste('OR[',1,',',1:length(data.ls$mk),']',sep='')]
	        data.nma.orsp<-data.nma[paste('OR[',2,',',1:length(data.ls$mk),']',sep='')]
	        data.nma.dor<-data.nma[paste('DOR[',1:length(data.ls$mk),']',sep='')]
            data.nma.sindex<-data.nma[paste('S[',1:length(data.ls$mk),']',sep='')]
            data.nma.sindex_se<-data.nma[paste('S_se[',1:length(data.ls$mk),']',sep='')]
            data.nma.sindex_sp<-data.nma[paste('S_sp[',1:length(data.ls$mk),']',sep='')]
	        mu.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.musp,mean))),
	                            ub_sp=array(unlist(lapply(data.nma.musp,func.ub))),
	                            lb_sp=array(unlist(lapply(data.nma.musp,func.lb))),
	                            mean_se=array(unlist(lapply(data.nma.muse,mean))),
	                            ub_se=array(unlist(lapply(data.nma.muse,func.ub))),
	                            lb_se=array(unlist(lapply(data.nma.muse,func.lb))))
	        rr.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.rrsp,mean))),
	                            ub_sp=array(unlist(lapply(data.nma.rrsp,func.ub))),
	                            lb_sp=array(unlist(lapply(data.nma.rrsp,func.lb))),
	                            mean_se=array(unlist(lapply(data.nma.rrse,mean))),
	                            ub_se=array(unlist(lapply(data.nma.rrse,func.ub))),
	                            lb_se=array(unlist(lapply(data.nma.rrse,func.lb))))
	        or.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.orsp,mean))),
	                            ub_sp=array(unlist(lapply(data.nma.orsp,func.ub))),
	                            lb_sp=array(unlist(lapply(data.nma.orsp,func.lb))),
	                            mean_se=array(unlist(lapply(data.nma.orse,mean))),
	                            ub_se=array(unlist(lapply(data.nma.orse,func.ub))),
	                            lb_se=array(unlist(lapply(data.nma.orse,func.lb))))
	        dor.info<-data.frame(mean=array(unlist(lapply(data.nma.dor,mean))),
	                             ub=array(unlist(lapply(data.nma.dor,func.ub))),
	                             lb=array(unlist(lapply(data.nma.dor,func.lb))))
            sindex.info<-data.frame(mean=array(unlist(lapply(data.nma.sindex,mean))),
	                                ub=array(unlist(lapply(data.nma.sindex,func.ub))),
	                                lb=array(unlist(lapply(data.nma.sindex,func.lb))))
            sindex_se.info<-data.frame(mean=array(unlist(lapply(data.nma.sindex_se,mean))),
	                                   ub=array(unlist(lapply(data.nma.sindex_se,func.ub))),
	                                   lb=array(unlist(lapply(data.nma.sindex_se,func.lb))))
            sindex_sp.info<-data.frame(mean=array(unlist(lapply(data.nma.sindex_sp,mean))),
	                                   ub=array(unlist(lapply(data.nma.sindex_sp,func.ub))),
	                                   lb=array(unlist(lapply(data.nma.sindex_sp,func.lb))))
	        results.nma<-list(marker=mk.names[mk.names %in% data.ls$mk],
	                          mu=mu.info,
	                          rr=rr.info,
	                          or=or.info,
	                          dor=dor.info,
                              sindex=sindex.info,
                              sindex_se=sindex_se.info,
                              sindex_sp=sindex_sp.info)
	   return(results.nma)
}

fixlen.str<-function(string.arr=NULL){
                                     	 max.len<-max(nchar(as.character(string.arr)))
	                                     rule<-paste('% -',as.character(max.len),'s',sep='')
	                                     string.arr.fix<-sprintf(rule,string.arr)
	                                     return(string.arr.fix)
}
fit.extract<-function(fit.data=NULL,filename.save=NULL){
	                                                       results.sum<-summary(fit.data)[[1]]
	                                                       results.simdata<-fit.data@sim$samples[[1]]
	                                                       save(results.sum,results.simdata,file=filename.save)
	                                                       return(results.simdata)
}

stan.code.nma<-'
data{
     int N;  //number of comparison??? - 121
     int Nt; //number of test - 10
     int Ns; //number of study - 72
     int TP[N];
     int Dis[N];  //diseased
     int TN[N];
     int NDis[N]; //non-diseased
     int Study[N];
     int Test[N];
}
parameters{
           matrix[2, Nt] logitmu;
           vector[Ns] nu[2];
           matrix[Ns, Nt] delta[2];
           vector<lower=0>[Nt] tau[2]; //*
           vector<lower=0>[2] sigmab; 
           real<lower=-1, upper=1> rho;
}
transformed parameters{
                       matrix[Ns, 2] p_i[Nt];
                       matrix[2, Nt] MU;
                       matrix[2, Nt] RR;
                       matrix[2, Nt] OR;
                       vector[Nt] DOR;
                       vector[Nt] S;
                       vector[Nt] S_se;
                       vector[Nt] S_sp;
                       matrix[Nt, Nt] A;
                       matrix[Nt, Nt] B;
                       matrix[Nt, Nt] C;
                       matrix[Nt, Nt] A_se;
                       matrix[Nt, Nt] B_se;
                       matrix[Nt, Nt] C_se;
                       matrix[Nt, Nt] A_sp;
                       matrix[Nt, Nt] B_sp;
                       matrix[Nt, Nt] C_sp;
                       vector<lower=0>[Nt] tausq[2];
                       vector<lower=0>[2] sigmabsq;
                       matrix[Nt, Nt] sigmasq[2];
                       matrix[Nt, Nt] rhow[2];

    for (i in 1:Ns){
        for (j in 1:2){
            for (k in 1:Nt)
                p_i[k][i,j] = inv_logit(logitmu[j,k] +  nu[j][i] + delta[j][i,k]);
        }
    }
 
    for (j in 1:2){
        for (k in 1:Nt){
                         MU[j,k] = mean(col(p_i[k], j));
        }
        tausq[j] = (tau[j]).*(tau[j]);
    }

    for (j in 1:2){
        for (k in 1:Nt){
                         RR[j, k] = MU[j, k]/MU[j, 1]; 
                         OR[j, k] = (MU[j, k]*(1 - MU[j, 1]))/(MU[j, 1]*(1 - MU[j, k]));
         }
    }

    for (l in 1:Nt){
                     DOR[l] = (MU[1, l]*MU[2, l])/((1 - MU[1, l])*(1 - MU[2, l]));

        for(m in 1:Nt){
                        A[l, m] = if_else((MU[1, l] > MU[1, m]) && (MU[2, l] > MU[2, m]), 1, 0);
                        B[l, m] = if_else((MU[1, l] < MU[1, m]) && (MU[2, l] < MU[2, m]), 1, 0);
                        C[l, m] = if_else((MU[1, l] == MU[1, m]) && (MU[2, l] == MU[2, m]), 1, 0);

                        A_se[l, m] = if_else((MU[1, l] > MU[1, m]), 1, 0);
                        B_se[l, m] = if_else((MU[1, l] < MU[1, m]), 1, 0);
                        C_se[l, m] = if_else((MU[1, l] == MU[1, m]), 1, 0);

                        A_sp[l, m] = if_else((MU[2, l] > MU[2, m]), 1, 0);
                        B_sp[l, m] = if_else((MU[2, l] < MU[2, m]), 1, 0);
                        C_sp[l, m] = if_else((MU[2, l] == MU[2, m]), 1, 0);
        }

        S[l] = (2*sum(row(A, l)) + sum(row(C, l)))/(2*sum(row(B, l)) + sum(row(C, l)));
        S_se[l] = (2*sum(row(A_se, l)) + sum(row(C_se, l)))/(2*sum(row(B_se, l)) + sum(row(C_se, l)));
        S_sp[l] = (2*sum(row(A_sp, l)) + sum(row(C_sp, l)))/(2*sum(row(B_sp, l)) + sum(row(C_sp, l)));
    }
    
    sigmabsq = (sigmab).*(sigmab);

    for (j in 1:2){
        for (k in 1:Nt){
            for (l in 1:Nt){
                             sigmasq[j][k,l] = (sigmabsq[j] + tausq[j][k])*((sigmabsq[j] + tausq[j][l]));
                             rhow[j][k,l] = sigmabsq[j]/sqrt(sigmasq[j][k,l]);
            }
        }
    }

}
model{
	   //Priors
       for (j in 1:2){
                       logitmu[j] ~ normal(0, 5);
		               tau[j] ~ cauchy(0, 2.5);
       }

         sigmab ~ cauchy(0, 2.5);
	     rho ~ uniform(-1, 1);
         nu[2] ~ normal(0, sigmab[2]);
      
       for (i in 1:Ns){
                       nu[1][i] ~ normal((sigmab[1]/sigmab[2])*rho*nu[2][i], sqrt(sigmabsq[1]*(1 - (rho*rho))));
          for (j in 1:2){
              for (k in 1:Nt)
                              delta[j][i,k] ~ normal(0, tau[j][k]);
        }
    }

    for (n in 1:N){
                    TP[n] ~ binomial(Dis[n], p_i[Test[n]][Study[n], 1]);
                    TN[n] ~ binomial(NDis[n], p_i[Test[n]][Study[n], 2]);
    }

}
generated quantities{
    
    vector[2*N] loglik;

    for (n in 1:N)
                   loglik[n] = binomial_lpmf(TN[n]| NDis[n], p_i[Test[n]][Study[n], 1]);

    for (n in (N+1):(2*N))
                   loglik[n] = binomial_lpmf(TN[n-N]| NDis[n-N], p_i[Test[n-N]][Study[n-N], 2]);

}
'

mk.names=c('qSOFA',
           'PCT',
           'presepsin',
           'CRP',
           'CD64',
           'IL-6',
           'sTREM-1',
           'LBP',
           'SIRS',
           'SOFA')

filename<-'Sepsis_nma_data.csv'

para.ls<-getpara(filename=filename,names.mk=mk.names)
para<-para.ls$para
data<-para.ls$data

para.as.ls.nma<-list(N=nrow(para),
                     Nt=length(unique(para$Test)),
                     Ns=length(unique(para$Study)),
                     TP=para$TP,
                     Dis=para$Dis,
                     TN=para$TN,
                     NDis=para$NDis,
                     Study=as.numeric(para$Study),
                     Test=para$Test)

stan.model.nma<-stan(model_code=stan.code.nma,
                     data=para.as.ls.nma, 
                     iter=1, 
                     warmup=0, 
                     chains=2)

start_time <- Sys.time()
fit.nma<-stan(fit=stan.model.nma,
              data=para.as.ls.nma,
              iter=20000,
              chains=4,
              cores=4,
              control=list(adapt_delta=0.99,
                           stepsize=0.01,
                           max_treedepth=17))
end_time <- Sys.time()
end_time - start_time

filename.simdata.save<-'fit_nma_simdata_qSOFAasreference.RData'
filename.results.save<-'out_nma_sepsis3_qSOFAasreference.csv'
results.nma<-fit.extract(fit.data=fit.nma,filename.save=filename.simdata.save)
results.nma.95<-get.nma(data.ls=list(mk=mk.names[unique(para$Test)],
                                         data.nma=results.nma),mk.names=mk.names,
                            conf=95)
write.csv(results.nma.95,filename.results.save)

rownames(results.nma.95$mu)=c('qSOFA','PCT','presepsin','CRP','CD64','IL-6','sTREM-1','LBP','SIRS','SOFA')
results.nma.95$mu

rownames(results.nma.95$dor)=c('qSOFA','PCT','presepsin','CRP','CD64','IL-6','sTREM-1','LBP','SIRS','SOFA')
results.nma.95$dor

rownames(results.nma.95$sindex)=c('qSOFA','PCT','presepsin','CRP','CD64','IL-6','sTREM-1','LBP','SIRS','SOFA')
results.nma.95$sindex

rownames(results.nma.95$sindex_se)=c('qSOFA','PCT','presepsin','CRP','CD64','IL-6','sTREM-1','LBP','SIRS','SOFA')
results.nma.95$sindex_se

rownames(results.nma.95$sindex_sp)=c('qSOFA','PCT','presepsin','CRP','CD64','IL-6','sTREM-1','LBP','SIRS','SOFA')
results.nma.95$sindex_sp

get.nma.1mk<-function(data.ls=NULL,mk.names=NULL,n.study=1,conf=95){
	           data.nma<-data.ls$data.nma
	           func.ub<-NULL
	           func.lb<-NULL
	           if(conf==95){
		                      func.ub<-get.ub.95
		                      func.lb<-get.lb.95
	           } else if(conf==90){
		                             func.ub<-get.ub.90
		                             func.lb<-get.lb.90
	           }
	           data.nma.muse<-data.nma[paste('MU[',1,',',1:length(data.ls$mk),']',sep='')]
	           data.nma.musp<-data.nma[paste('MU[',2,',',1:length(data.ls$mk),']',sep='')]
	           data.nma.rrse<-data.nma[paste('RR[',1,',',1:length(data.ls$mk),']',sep='')]
	           data.nma.rrsp<-data.nma[paste('RR[',2,',',1:length(data.ls$mk),']',sep='')]
	           data.nma.orse<-data.nma[paste('OR[',1,',',1:length(data.ls$mk),']',sep='')]
	           data.nma.orsp<-data.nma[paste('OR[',2,',',1:length(data.ls$mk),']',sep='')]
	           data.nma.dor<-data.nma[paste('DOR[',1:length(data.ls$mk),']',sep='')]
	           data.nma.sindex<-data.nma[paste('S[',1:length(data.ls$mk),']',sep='')]
	           data.nma.sindex_se<-data.nma[paste('S_se[',1:length(data.ls$mk),']',sep='')]
	           data.nma.sindex_sp<-data.nma[paste('S_sp[',1:length(data.ls$mk),']',sep='')]
	           data.nma.beta1se<-data.nma[paste('beta1[1,',1:length(data.ls$mk),']',sep='')]
	           data.nma.beta1sp<-data.nma[paste('beta1[2,',1:length(data.ls$mk),']',sep='')]
	           mu.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.musp,mean))),
	                               ub_sp=array(unlist(lapply(data.nma.musp,func.ub))),
	                               lb_sp=array(unlist(lapply(data.nma.musp,func.lb))),
	                               mean_se=array(unlist(lapply(data.nma.muse,mean))),
	                               ub_se=array(unlist(lapply(data.nma.muse,func.ub))),
	                               lb_se=array(unlist(lapply(data.nma.muse,func.lb))))
	           rr.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.rrsp,mean))),
	                               ub_sp=array(unlist(lapply(data.nma.rrsp,func.ub))),
	                               lb_sp=array(unlist(lapply(data.nma.rrsp,func.lb))),
	                               mean_se=array(unlist(lapply(data.nma.rrse,mean))),
	                               ub_se=array(unlist(lapply(data.nma.rrse,func.ub))),
	                               lb_se=array(unlist(lapply(data.nma.rrse,func.lb))))
	           or.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.orsp,mean))),
	                               ub_sp=array(unlist(lapply(data.nma.orsp,func.ub))),
	                               lb_sp=array(unlist(lapply(data.nma.orsp,func.lb))),
	                               mean_se=array(unlist(lapply(data.nma.orse,mean))),
	                               ub_se=array(unlist(lapply(data.nma.orse,func.ub))),
	                               lb_se=array(unlist(lapply(data.nma.orse,func.lb))))
	           dor.info<-data.frame(mean=array(unlist(lapply(data.nma.dor,mean))),
	                                ub=array(unlist(lapply(data.nma.dor,func.ub))),
	                                lb=array(unlist(lapply(data.nma.dor,func.lb))))
	           sindex.info<-data.frame(mean=array(unlist(lapply(data.nma.sindex,mean))),
	                                   ub=array(unlist(lapply(data.nma.sindex,func.ub))),
	                                   lb=array(unlist(lapply(data.nma.sindex,func.lb))))
	           sindex_se.info<-data.frame(mean=array(unlist(lapply(data.nma.sindex_se,mean))),
	                                      ub=array(unlist(lapply(data.nma.sindex_se,func.ub))),
	                                      lb=array(unlist(lapply(data.nma.sindex_se,func.lb))))
	           sindex_sp.info<-data.frame(mean=array(unlist(lapply(data.nma.sindex_sp,mean))),
	                                      ub=array(unlist(lapply(data.nma.sindex_sp,func.ub))),
	                                      lb=array(unlist(lapply(data.nma.sindex_sp,func.lb))))
	           beta1.info<-data.frame(mean_sp=array(unlist(lapply(data.nma.beta1sp,mean))),
	                                  ub_sp=array(unlist(lapply(data.nma.beta1sp,func.ub))),
	                                  lb_sp=array(unlist(lapply(data.nma.beta1sp,func.lb))),
	                                  mean_se=array(unlist(lapply(data.nma.beta1se,mean))),
	                                  ub_se=array(unlist(lapply(data.nma.beta1se,func.ub))),
	                                  lb_se=array(unlist(lapply(data.nma.beta1se,func.lb))))
	           results.nma<-list(marker=mk.names[mk.names %in% data.ls$mk],
	                             mu=mu.info,
	                             rr=rr.info,
	                             or=or.info,
	                             dor=dor.info,
	                             sindex=sindex.info,
	                             sindex_se=sindex_se.info,
	                             sindex_sp=sindex_sp.info,
	                             beta1=beta1.info)
	return(results.nma)
}

stan.code.1mk<-'
                data{
                      int N;         //number of comparison (n=121)
                      int Nt;        //number of test (n=7)
                      int Ns;        //number of study (n=107)
                      int TP[N];
                      int Dis[N];    //diseased
                      int TN[N];
                      int NDis[N];   //non-diseased
                      int Study[N];
                      int Test[N];
                      int Covar1[N]; //one additional variable
                    }

                parameters{
                            matrix[2,Nt] logitmu;
                            vector[Ns] nu[2];
                            matrix[Ns,Nt] delta[2];
                            vector<lower=0>[2] tau;
                            vector<lower=0>[2] sigmab; 
                            real<lower=-1,upper=1> rho;
                            matrix[2,Nt] beta1;
                        }

                transformed parameters{
                                        matrix[Ns,2] p_i[Nt];
                                        matrix[2,Nt] MU;
                                        matrix[2,Nt] RR;
                                        matrix[2,Nt] OR;
                                        vector[Nt] DOR;
                                        vector[Nt] S;
                                        vector[Nt] S_se;
                                        vector[Nt] S_sp;
                                        matrix[Nt,Nt] A;
                                        matrix[Nt,Nt] B;
                                        matrix[Nt,Nt] C;
                                        matrix[Nt, Nt] A_se;
                                        matrix[Nt, Nt] B_se;
                                        matrix[Nt, Nt] C_se;
                                        matrix[Nt, Nt] A_sp;
                                        matrix[Nt, Nt] B_sp;
                                        matrix[Nt, Nt] C_sp;

                                        vector<lower=0>[2] tausq;
                                        vector<lower=0>[2] sigmabsq;

                                        matrix[Nt,Nt] sigmasq[2];
                                        matrix[Nt,Nt] rhow[2];

                for (i in 1:Ns){
                     for (j in 1:2){
                          for (k in 1:Nt){
                                           p_i[k][i,j]=inv_logit(logitmu[j,k]+nu[j][i]+delta[j][i,k]+beta1[j,k]*Covar1[i]);
                           }
                      }
                 }

                for (j in 1:2){
                     for (k in 1:Nt){
                                      MU[j,k]=mean(col(p_i[k],j));
                      }
                }
                  
                for (j in 1:2){
                     for (k in 1:Nt){
                                      RR[j,k]=MU[j,k]/MU[j,1]; 
                                      OR[j,k]=(MU[j,k]*(1-MU[j,1]))/(MU[j,1]*(1-MU[j,k]));
                       }
                 }

                for (l in 1:Nt){
                                 DOR[l]=(MU[1,l]*MU[2,l])/((1-MU[1,l])*(1-MU[2,l]));
                     for(m in 1:Nt){
                                     A[l,m]=if_else((MU[1,l]>MU[1,m]) && (MU[2,l]>MU[2,m]),1,0);
                                     B[l,m]=if_else((MU[1,l]<MU[1,m]) && (MU[2,l]<MU[2,m]),1,0);
                                     C[l,m]=if_else((MU[1,l]==MU[1,m]) && (MU[2,l]==MU[2,m]),1,0);
                                     
                                     A_se[l, m] = if_else((MU[1, l] > MU[1, m]), 1, 0);
                                     B_se[l, m] = if_else((MU[1, l] < MU[1, m]), 1, 0);
                                     C_se[l, m] = if_else((MU[1, l] == MU[1, m]), 1, 0);
 
                                     A_sp[l, m] = if_else((MU[2, l] > MU[2, m]), 1, 0);
                                     B_sp[l, m] = if_else((MU[2, l] < MU[2, m]), 1, 0);
                                     C_sp[l, m] = if_else((MU[2, l] == MU[2, m]), 1, 0);
                      }
                        S[l]=(2*sum(row(A,l))+sum(row(C,l)))/(2*sum(row(B,l))+sum(row(C,l)));
                        S_se[l] = (2*sum(row(A_se, l)) + sum(row(C_se, l)))/(2*sum(row(B_se, l)) + sum(row(C_se, l)));
                        S_sp[l] = (2*sum(row(A_sp, l)) + sum(row(C_sp, l)))/(2*sum(row(B_sp, l)) + sum(row(C_sp, l)));
                }
                 tausq=(tau).*(tau);  
                 sigmabsq=(sigmab).*(sigmab);

                for (j in 1:2){
                     for (k in 1:Nt){
                          for (l in 1:Nt){
                                           sigmasq[j][k,l]=(sigmabsq[j]+tausq[j])*((sigmabsq[j]+tausq[j]));
                                           rhow[j][k,l]=sigmabsq[j]/sqrt(sigmasq[j][k,l]);
                            }
                       }
                 }
            }

                model{
                      //Priors
                      for (j in 1:2){
                                      logitmu[j]~normal(0,5);
                                      tau[j]~cauchy(0,2.5);
                                      sigmab[j]~gamma(5,1);
                           for (k in 1:Nt){
                                            beta1[j,k]~normal(0,5);
                           }            
                      }
                       rho~uniform(-1, 1);
                       nu[2]~normal(0,sigmab[2]);

                      for (i in 1:Ns){
                                       nu[1][i]~normal((sigmab[1]/sigmab[2])*rho*nu[2][i],sqrt(sigmabsq[1]*(1-(rho*rho))));
                           for (j in 1:2){
                               for (k in 1:Nt){
                                                delta[j][i,k]~normal(0,tau[j][k]);
                              }
                         }
                      }

                     for (n in 1:N){
                                     TP[n]~binomial(Dis[n],p_i[Test[n]][Study[n],1]);
                                     TN[n]~binomial(NDis[n],p_i[Test[n]][Study[n],2]);
                      }

                    }

                generated quantities{
                                      vector[2*N] loglik;

                                      for (n in 1:N){
                                                      loglik[n]=binomial_lpmf(TN[n]|NDis[n],p_i[Test[n]][Study[n],1]);
                                       }
                                      for (n in (N+1):(2*N)){
                                                              loglik[n]=binomial_lpmf(TN[n-N]|NDis[n-N],p_i[Test[n-N]][Study[n-N],2]);
                                       }
                }
'

para.as.ls.1mk<-list(N=nrow(para),
                     Nt=length(unique(para$Test)),
                     Ns=length(unique(para$Study)),
                     TP=para$TP,
                     Dis=para$Dis,
                     TN=para$TN,
                     NDis=para$NDis,
                     Study=as.numeric(para$Study),
                     Test=para$Test,
                     Covar1=data$is.prevalence05)

stan.model.1mk<-stan(model_code=stan.code.1mk, 
                     data=para.as.ls.1mk, 
                     iter=1, 
                     warmup=0, 
                     chains=1)

start_time <- Sys.time()
fit.1mk<-stan(fit=stan.model.1mk,
              data=para.as.ls.1mk,
              iter=20000,
              chains=4,
              cores=4,
              control=list(adapt_delta=0.99,
                           stepsize=0.01,
                           max_treedepth=17))
end_time <- Sys.time()
end_time - start_time

filename.simdata.save<-'fit_1mk_simdata_sepsis_isprevalence05.RData'
filename.results.save<-'out_1mk_sepsis_isprevalence05.csv'
results.1mk<-fit.extract(fit.data=fit.1mk,filename.save=filename.simdata.save)
results.1mk.95<-get.nma.1mk(data.ls=list(mk=mk.names[unique(para$Test)],
                                         data.nma=results.1mk),mk.names=mk.names,
                            conf=95)
write.csv(results.1mk.95,filename.results.save)

rownames(results.1mk.95$mu)=c('qSOFA', 'PCT', 'presepsin', 'CRP', 'CD64', 'IL-6', 'sTREM-1', 'LBP', 'SIRS', 'SOFA')
results.lmk.95$mu

rownames(results.1mk.95$dor)=c('qSOFA', 'PCT', 'presepsin', 'CRP', 'CD64', 'IL-6', 'sTREM-1', 'LBP', 'SIRS', 'SOFA')
results.1mk.95$dor
