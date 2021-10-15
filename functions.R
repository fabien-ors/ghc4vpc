#
#                          GHC4VPC
#
#          Geostatistical Hierarchical Clustering
#                        applied to
#             Vertical Sand Proportion Curves
#
# Version 1.1
# 2021, October 15th
#
# Copyright (c) 2021 - Thomas Romary & Fabien Ors
#
# LICENSE: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#
# REFERENCES:
#
# Romary T., Ors F., Rivoirard J., Deraisme J. 
# Unsupervised classification of multivariate geostatistical data: Two algorithms.
# Comput. Geosci. 85, 96–103 (2015).
# https://doi.org/10.1016/j.cageo.2015.05.019
#
# Bubnova, A., Ors, F., Rivoirard, J. et al.
# Automatic Determination of Sedimentary Units from Well Data.
# Math Geosci 52, 213–231 (2020).
# https://doi.org/10.1007/s11004-019-09793-w
#

#################################################################################
# Création d'un voisinage type "VPC"

create_vois=function(data) {
  n=dim(data)[1]
  m=matrix(0,n,n)
  for(i in 1:n-1) {
    m[i,i+1]=1
  }
  return(m)
} 

####################################################################################
# Calcule la matrice de distance euclidienne au carré pour les données de data
# data : data.frame de données
# w : vecteur de poids (size = nb variables)
# cont : vecteur de liste des variables continues
# cat : vecteur de liste des variables categorielles
# Pour l'instant, distance euclidienne au carré mais ajout
# facile d'autres distances

dis = function(data,w,cont,cat) {
  m=dim(data)[1]
  dis1=matrix(0,m,m)
  maxdis=0
  for (i in 2:m) {
    k=1
    for (j in cont){
      dis1[i,1:(i-1 )]=dis1[i,1:(i-1)]+w[k]*(data[i,j]-data[1:(i-1),j])^2
      k=k+1
    }
    for (j in cat) {
      dis1[i,1:(i-1)]=dis1[i,1:(i-1)]+w[k]*as.numeric(data[i,j]!=data[1:(i-1),j])
    }
  }
  # upper triangle
  return(t(dis1))
}

############################################################
# Clustering Hiérarchique Géostatistique
# dis1 : matrice de distance
# vois : matrice de voisinage
# method : critère de linkage parmi :
#          'single','complete','average','weighted_average',
#          'ward','wardp'
# keep : nombre de niveaux conservés (50 par défaut)

geohier = function(dis1,vois,method,keep = 50) {
  k=1
  m=dim(dis1)[1]
  clus=matrix(seq(1:m))
  nbclus = rep(1,m)
  #rq on ne conserve que les keep derniers niveaux
  clus_final=NULL
  heights = NULL # hauteurs (en vue dendrogramme)
  merge = NULL  
  # Paramètres de linkage selon formule de Lance-Williams
  link_param = switch(method,
                      single = c(0.5,0.5,0,-0.5),
                      complete = c(0.5,0.5,0,0.5),
                      average = c(1,1,0,0), 
                      weighted_average = c(0.5,0.5,-0.25,0),
                      ward = c(1,1,1,0),
                      wardp = c(1,1,1,0))
  
  while (k!=m-1) {
    dis2=dis1
    dis2[vois!=1]=10000
    mindis=min(dis2)
    tata=which(dis2==mindis)[1]
    toto=c(tata%%(m-k+1),tata%/%(m-k+1)+1)
    mc=toto[1]+(m-k+1)*(toto[1]==0) 
    ml=toto[2]-(toto[1]==0)         
    uclus=unique(clus)
    vecpos=as.numeric(vois[mc,]|vois[ml,]|vois[,mc]|vois[,ml])[-c(mc,ml)]
    dc = c(dis1[1:mc,mc],dis1[mc,(mc+1):(m-k+1)])
    if (ml == (m-k+1)) {
      dl = dis1[,ml]
    }
    else {
      dl = c(dis1[1:ml,ml],dis1[ml,(ml+1):(m-k+1)])
    }
    if (method == 'average') { 
      nc = sum(clus == uclus[mc])
      nl = sum(clus == uclus[ml])
      ns = nc+nl
      link_param = c(nc/ns,nl/ns,0,0)
    }
    if (method != 'ward' && method != 'wardp') {
      dispos=(link_param[1]*dc+
              link_param[2]*dl+
              link_param[3]*mindis+
              link_param[4]*abs(dc-dl))[-c(mc,ml)] # L-W formula
    }
    if (method == 'ward' || method == 'wardp') {
      ns = nbclus[mc]+nbclus[ml]+nbclus           # ns = ni+nj+nk
      # link param pour ward+ : vecteurs !
      link_param = cbind((nbclus[mc]+  nbclus)/ns, # ai = (ni+nk)/(ns)
                         (nbclus[ml]+  nbclus)/ns, # aj = (nj+nk)/(ns)
                         ifelse(method == 'ward',
                                -nbclus/ns,       # b  = -nk/(ns) (ward standard)
                                 nbclus/ns),      # b  = +nk/(ns) (ward plus)
                         rep(0,length(nbclus))) 
      dispos=(link_param[,1]*dc+
              link_param[,2]*dl+
              link_param[,3]*mindis+
              link_param[,4]*abs(dc - dl))[-c(mc,ml)] # L-W formula
      nbclus[mc] = (nbclus[mc] + nbclus[ml])
      nbclus = nbclus[-ml]
    }
    clus[clus==uclus[ml]]=uclus[mc]
    if (k>m-keep-1) {
      clus_final = cbind(clus_final,clus)
      heights = c(heights,mindis)
      merge = rbind(merge,c(mc,ml))
    }
    vois=vois[-ml,]
    vois=vois[,-ml]
    dis1=dis1[-ml,]
    dis1=dis1[,-ml]
    if (mc==1) {
      vois[mc,2:(m-k)]=vecpos
      dis1[mc,2:(m-k)]=dispos
    }
    else if (mc==(m-k)) {
      vois[1:(mc-1),mc]=vecpos
      dis1[1:(mc-1),mc]=dispos
    }
    else {
      vois[1:(mc-1),mc]=vecpos[1:(mc-1)]
      vois[mc,(mc+1):(m-k)]=vecpos[mc:(m-k-1)]
      dis1[1:(mc-1),mc]=dispos[1:(mc-1)]
      dis1[mc,(mc+1):(m-k)]=dispos[mc:(m-k-1)]
    }
    k=k+1
    print(k)
  }
  mc=1
  ml=2
  merge = rbind(merge,c(mc,ml))
  uclus=unique(clus)
  heights = c(heights,dis1[mc,ml])
  clus[clus==uclus[ml]]=uclus[mc]
  clus_final=cbind(clus_final,clus)
  print(k+1)  
  return(list(clust = clus_final, heights = heights, merge = merge, method = method))
}

##########################################################
# Post traitement / Affichage
# transforme les numéros de clusters
# en nombre entre 1 et nbclus
# data : tableau de données
# clus_final : objet clust
# nbclus : nombre de clust

postclus = function(clus_final,nbclus,keep = 50) {
  tri=sort(unique(clus_final$clust[,keep+1-nbclus]),index=T)
  color=clus_final$clust[,keep+1-nbclus]
  for (i in tri$ix) {
    color[color==tri$x[i]]=tri$ix[i]
  }
  return(color)
}

write.clust = function(clus,clusname) {
  system(paste('mkdir ',clusname,sep=''))
  write.csv(clus$clust,paste(clusname,'/clus.csv',sep=''))
  write.csv(clus$heights,paste(clusname,'/heights.csv',sep=''))
  write.csv(clus$merge,paste(clusname,'/merge.csv',sep=''))
  write.csv(clus$method,paste(clusname,'/method.csv',sep=''))
  return()
}

read.clust = function(clusname) {
  return(list(clust = read.csv(paste(clusname,'/clus.csv',sep=''))[,-1],
              heights = read.csv(paste(clusname,'/heights.csv',sep=''))[,-1],
              merge = read.csv(paste(clusname,'/merge.csv',sep=''))[,-1],
              method = read.csv(paste(clusname,'/method.csv',sep=''))[,-1]))
}

