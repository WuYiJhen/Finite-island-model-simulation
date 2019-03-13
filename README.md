# Finite-island-model-simulation




#' Simulate stochastic process under IAM-FIM
#' IAM.FIM.Simu(x,m,u,g) is a function of obtaining the MLE of Shannon entropy and heterozygosity of total population and subpopulation in each generation under IAM-FIM.
#' @param x is the matrix of species of the first generation. Row length is the number of diploid individuals in an idealized population. Column length is the number of idealized subpopulations.
#' @param m is the dispersal (or migration) rate per generation.
#' @param u is the mutation rate per generation.
#' @param g is the generations of the stochastic process.
#' @return the MLE of Shannon entropy and heterozygosity of total population and subpopulation in each generation.

IAM.FIM.Simu <- function(x,m,u,g){
  cppFunction('NumericMatrix cppCount(NumericVector name, NumericMatrix x){
    int n=name.size();
    int nrow=x.nrow();
    int ncol=x.ncol();
    NumericMatrix res(n,ncol);
    for(int k=0;k<ncol;++k){
      for(int i=0;i<n;++i){
        int a=0;
        for(int j=0;j<nrow;++j){
          if(name[i]==x(j,k)){
            a+=1;
          }
        }
        res(i,k)=a;
      }
    }
    return res;
  }')
  cppFunction('NumericVector hmle(NumericMatrix x){
    int s=x.nrow();
    int n=x.ncol();
    int z=0;
    for(int i=0;i<s;++i){
      z+=x(i,0);
    }
    NumericMatrix X(s,n);
    NumericMatrix X1(s,n);
    for(int j=0;j<n;++j){
      for(int i=0;i<s;++i){
        X(i,j)=(double)x(i,j)/z;
        X1(i,j)=x(i,j)*x(i,j);
      }
    }
    NumericVector y(s);
    NumericVector y1(s);
    for(int i=0;i<s;++i){
      for(int j=0;j<n;++j){
        y[i]+=(double)X(i,j)/n;
        y1[i]+=X1(i,j);
      }
    }
    double ht=0;
    double ht2=0;
    for(int i=0;i<s;++i){
      ht+=(double)y[i]*log(y[i]);
      ht2+=(double)pow(y[i],2);
    }
    NumericVector hj(n);
    for(int j=0;j<n;++j){
      for(int i=0;i<s;++i){
        if(X(i,j)>0){
          hj[j]+=(double)X(i,j)*log(X(i,j))/n;
        }
      }
    }
    double hs=(double)std::accumulate(hj.begin(),hj.end(),0.0);
    double hs2=(double)std::accumulate(y1.begin(),y1.end(),0.0);
    NumericVector res(4);
    int kk=n*z*z;
    res[0]=-ht;
    res[1]=-hs;
    res[2]=1-ht2;
    res[3]=1-(double)hs2/kk;
    return res;
  }')
  cppFunction('NumericMatrix cppIAM_FIM(NumericMatrix x, double m ,double u){
    int ncol=x.ncol();
    int nrow=x.nrow();
    NumericMatrix mutate(nrow,ncol);
    NumericMatrix migrate(nrow,ncol);
    for(int i=0;i<ncol;++i){
      mutate(_,i)=rbinom(nrow,1,u);
      migrate(_,i)=rbinom(nrow,1,m);
    }
    int ng=max(x)+1;
    for(int j=0;j<ncol;++j){
      for(int i=0;i<nrow;++i){
        if(mutate(i,j)>0){
          x(i,j)=ng;
          ng+=1;
        }
      }
    }
    NumericMatrix Ng(nrow,ncol);
    NumericVector il(ncol);
    for(int j=0;j<ncol;++j){
      for(int i=0;i<nrow;++i){
        int k=rand()%ncol;
        if(migrate(i,j)>0){
          while(k==j){
            k=rand()%ncol;
          }
          Ng(il[k],k)=x(i,j);
          il[k]+=1;
        }
        if(migrate(i,j)<1 && mutate(i,j)>0){
          Ng(il[j],j)=x(i,j);
          il[j]+=1;
        }
      }
    }
    for(int j=0;j<ncol;++j){
      NumericVector gg(nrow);
      int left=nrow-il[j];
      int naive=0;
      int a=0;
      for(int i=0;i<nrow;++i){
        if(mutate(i,j)==0 && migrate(i,j)==0){
          gg[a]=x(i,j);
          a+=1;
          naive+=1;
        }
      }
      for(int i=0;i<left;++i){
        int k1=rand()%naive;
        int index=il[j]+i;
        Ng(index,j)=gg[k1];
      }
    }
    return Ng;
  }')
  N <- nrow(x)
  n <- ncol(x)
  Hts <- matrix(0,4,(g+1))
  name <- unique(c(x))
  X <- cppCount(name,x)
  Hts[,1] <- hmle(X)
  for(i in 2:(g+1)){
    x1<- cppIAM_FIM(x,m,u)
    name <- unique(c(x1))
    X <- cppCount(name,x1)
    Hts[,i] <- hmle(X)
    x <- x1
  }
  return(list("Ht.1"=Hts[1,],"Hs.1"=Hts[2,],"Ht.2"=Hts[3,],"Hs.2"=Hts[4,]))
  } 
