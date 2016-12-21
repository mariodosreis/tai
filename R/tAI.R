##########################################################################
# This file contains a series of R functions for calculating the tRNA
# adaptation index (tAI), the S value (or correlation between tAI and the
# adjusted Nc) and other auxiliary functions.
#
# Copyright (c) 2003 Mario dos Reis
#
# This file forms part of codonR, a package for the analysis of codon
# usage in DNA coding sequences.
#
# codonR is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# Mario dos Reis
# m.reis@mail.cryst.bbk.ac.uk
# School of Crystallography
# Birkbeck College
# London WC1E 7HX
#
# version 0.2 (Nov 2003)
###########################################################################

##########################################################################
# Optimised Nc adjusted function:
##########################################################################

nc.adj <- function(nc, gc3) {

  a = -6.0    # a, b, and c have already been optimised
  b = 34.0
  c = 1.025
  
  x = a + gc3 + (b/(gc3^2 + (c - gc3)^2)) - nc
  return(x)
}

##########################################################################
# Function to calculate relative adaptiveness values
##########################################################################

# non-optimised s-values:
# s <- c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
# optimised s-values:
ops <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)

get.ws <- function(tRNA,  # tRNA gene copy number
                   s = ops,     # selective constraints
                   sking) # super kingdom: 0-eukaryota, 1-prokaryota
{

  p = 1 - s

  # initialise w vector
  W = NULL  # don't confuse w (lowercase) and W (uppercase)

  # obtain absolute adaptiveness values (Ws)
  for (i in seq(1, 61, by=4))
    W = c(W,
      p[1]*tRNA[i]   + p[5]*tRNA[i+1],     # INN -> NNT, NNC, NNA
      p[2]*tRNA[i+1] + p[6]*tRNA[i],       # GNN -> NNT, NNC
      p[3]*tRNA[i+2] + p[7]*tRNA[i],       # TNN -> NNA, NNG
      p[4]*tRNA[i+3] + p[8]*tRNA[i+2])     # CNN -> NNG

  # check methionine
  W[36] = p[4]*tRNA[36]
  
  # if bacteria, modify isoleucine ATA codon
  if(sking == 1) W[35] = p[9]

  # get rid of stop codons (11, 12, 15) and methionine (36)
  W = W[-c(11,12,15,36)]

  # get ws
  w = W/max(W)

  if(sum(w == 0) > 0) {
    ws <- w[w != 0] # zero-less ws
    gm <- exp(sum(log(ws))/length(ws)) # geometric mean
    w[w == 0] = gm # substitute 0-ws by gm
  }

  return(w)
}

#############################################################################
# After calculating all w's, the tRNA adaptation index (tAI) can then be
# calculated. x is a matrix with all the codon frequencies per ORF in any
# analysed genome
#############################################################################

get.tai <- function(x,w) {

  w = log(w)			#calculate log of w
  n = apply(x,1,'*',w)		#multiply each row of by the weights
  n = t(n)                      #transpose
  n = apply(n,1,sum)		#sum rows
  L = apply(x,1,sum)		#get each ORF length
  tAI = exp(n/L)		#get tai
  return(tAI)
}

############################################################################
# After calculating tAI, we can then calculate S, i.e. the correlation
# between tAI and g(GC3s) - Nc (aka nc.adj)
############################################################################

get.s <- function(tAI, nc, gc3) {

  nc.ad <- nc.adj(nc, gc3)
  ts <- cor(tAI, nc.ad, use='p')
  return(ts)
}

#############################################################################
# Statistical test for tAI
#############################################################################

ts.test <- function(m, ws, nc, gc3s, ts.obs, samp.size=500, n=1000) {

  # create a matrix of randomly permuted w-values:
  ws.permuted <- matrix(rep(ws, n), ncol = 60, byrow = T)
  ws.permuted <- t(apply(ws.permuted, 1, sample))

  # initialise 'translational selection' (ts) vector:
  ts = numeric(length(ws.permuted[,1]))
  # work with a smaller m matrix to reduce computation time:
  samp <- sample(nrow(m), samp.size)
  m.samp <- m[samp,]
  nc <- nc[samp]
  gc3s <- gc3s[samp]
  nc.ad <- nc.adj(nc, gc3s)

  # calculate simulated ts values:
  for(i in 1:length(ws.permuted[,1])) {
    tai <- get.tai(m.samp, ws.permuted[i,])
    co <- cor(nc.ad, tai, use = 'p')
    ts[i] <- co
  }

  # estimate the density function of ts values:
  ts.den <- density(ts)
  # bin width:
  ts.width <- ts.den$x[2] - ts.den$x[1]
  # p-values is:
  p.value <- 1 - sum(ts.width * ts.den$y[ts.den$x < ts.obs])
  # This method is approximate, it might happen that the above sum > 1
  # then, p.value needs to be corrected:
  if(p.value < 0) p.value <- 0

  return(list(p.value=p.value, ts.simulated=ts))
}

#############################################################################
# Nc plotting function
#############################################################################

nc.f <- function(x) -6 + x + (34/(x^2 + (1.025 - x)^2))

plot.nc <- function(nc, gc3s) {
  curve(nc.f, xlim = c(0, 1), ylim = c(20, 62),
        xlab = "GC3s", ylab = "Nc")
  points(nc ~ gc3s, pch = '+', cex = 0.5)
}
