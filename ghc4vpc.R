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

source("functions.R")

# Chargement d'une VPC (fichiers exporté par le logiciel
# Flumy : https://flumy.minesparis.psl.eu
# VPC Loranca (8 puits) : Stats1m
# VPC synthétique (20 puits / 3 unités) : 3u_20w
data=read.csv("3u_20w",header=T,sep=',',skip=25,nrows=85)
weight=c(1)

# Matrice des distances + voisinage
dis1=dis(data,weight,c(19),c())
vois=create_vois(data)

# Clustering Hiérarchique Géostatistique
clus=geohier(dis1,vois,'wardp')

# Affichage des hauteurs de clusters
plot(clus$heights)
X11()

# Affichage des clusters sur la VPC
nbclus=3
color=postclus(clus,nbclus)
rb=rainbow(nbclus)
z=data[[1]]
h=data[[19]]
ticks=unlist(c(0,lapply(seq(1,nbclus-1), function(x) { min(z[color>x]) }),max(z)+1))
ticks2=ticks*1.2
barplot(height=h,horiz=T,col=rb[color],xlim=c(0,1))
axis(side=2,at=ticks2,labels=as.character(ticks))

