library(corrplot)

#Functions to make the X design matrix from a matrix containing the lobe and hemisphere labels of the nodes

level_matrix<-function(num)
{
  levelmat=matrix(0,num,num)
  s=1
  for(i in 1:num)
  {
    for(j in 1:i)
    {
      levelmat[i,j]=s
      levelmat[j,i]=s
      s=s+1
    }
  }
  levelmat
}


make_X<-function(LobeHemiMatrix)
{
  V=nrow(LobeHemiMatrix)
  hem_names=sort(unique(LobeHemiMatrix$Hemisphere))
  lobe_names=sort(unique(LobeHemiMatrix$Lobe))
  N_hemi=length(hem_names)
  N_lobe=length(lobe_names)
  L=V*(V-1)/2
  
  h=level_matrix(N_hemi)
  l=level_matrix(N_lobe)
  
  X1=matrix(0,L,max(h))
  #col-names
  names1=NULL
  for(i in 1:N_hemi)
  {
    for(j in 1:i)
    {
      names1=c(names1,paste0("hemi",hem_names[i],".","hemi",hem_names[j]))
    }
  }
  colnames(X1)=names1
  X2=matrix(0,L,max(l))
  #col-names
  names2=NULL
  for(i in 1:N_lobe)
  {
    for(j in 1:i)
    {
      names2=c(names2,paste0("lobe",lobe_names[i],".","lobe",lobe_names[j]))
    }
  }
  colnames(X2)=names2
  
  i=0
  for(co in 2:V)
  {
    for(ro in 1:(co-1))
    {
      i=i+1
      
      hemi1=which(hem_names==LobeHemiMatrix$Hemisphere[co])
      hemi2=which(hem_names==LobeHemiMatrix$Hemisphere[ro])
      X1[i,h[hemi1,hemi2]]=1
      
      lobe1=which(lobe_names==LobeHemiMatrix$Lobe[co])
      lobe2=which(lobe_names==LobeHemiMatrix$Lobe[ro])
      X2[i,l[lobe1,lobe2]]=1
    }
  }
  
  X=cbind(X1,X2)[,-1]
  return(X)
}


#Function to convert adjacency matrix to a vector

matrix_2_vector<-function(a)  a[upper.tri(a)]


#Functions to plot the adjacency matrix

vector_2_matrix<-function(a,d=FALSE) 
{
  n=ceiling(sqrt(2*length(a)))-d
  C=matrix(0,nrow=n,ncol=n)
  C[upper.tri(C,diag = d)]=a
  C+t(C)-diag(diag(C))
}

show_connectome<-function(a)
{
  a=vector_2_matrix(a)
  corrplot::corrplot(a,is.corr=F)
}

