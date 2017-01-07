##########################################################################
# This file contains a series of R functions for calculating the tRNA
# adaptation index (tAI), the S value (or correlation between tAI and the
# adjusted Nc) and other auxiliary functions.
#
# Copyright (c) 2016-2003 Mario dos Reis
#
# This file forms part of tAI, a package for the analysis of codon
# usage in DNA coding sequences.
#
# tAI is free software; you can redistribute it and/or
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
# mariodosreis@gmail.com
#
# version 0.21 (Dec 2016)
###########################################################################

##########################################################################
# Adjusted Nc
##########################################################################
#' Adjusted effective number of codons (Nc)
#' @author Mario dos Reis
#' @export
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


#' Relative adaptiveness values
#' 
#' @author Mario dos Reis
#' 
#' @examples 
#' eco.ws <- get.ws(tRNA=ecolik12$trna, sking=1)
#' 
#' @export
get.ws <- function(tRNA,  # tRNA gene copy number
                   s = NULL,     # selective constraints
                   sking) # super kingdom: 0-eukaryota, 1-prokaryota
{
  # optimised s-values:
  if(is.null(s)) s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)

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
#' The tRNA adaptation index
#' 
#' Calculates the tRNA adaptation index (tAI) of dos Reis et al. (2003, 2004).
#' 
#' @param x an n by 60 matrix of relative codon frequencies for n open reading frames.
#' @param w a vector of length 61 of relative adaptiveness values for codons.
#' 
#' @details The tRNA adaptation index (tAI) is a measure of the level of 
#' co-adaptation between the set of tRNA genes and the codon usage bias of 
#' protein-coding genes in a given genome. STOP and methionine codons are 
#' ignored. The standard genetic code is assumed.
#' 
#' @return A vector of length n of tAI values. 
#' 
#' @author Mario dos Reis
#' 
#' @references 
#' dos Reis M., Wernisch L., and Savva R. (2003) Unexpected correlations 
#' between gene expression and codon usage bias from microarray data for the 
#' whole \emph{Escherichia coli} K-12 genome. \emph{Nucleic Acids Res.},
#' \bold{31:} 6976--85.
#' 
#' dos Reis M., Savva R., and Wernisch L. (2004) Solving the riddle of codon 
#' usage preferences: a test for translational selection. \emph{Nucleic Acids Res.},
#' \bold{32:} 5036--44.
#' 
#' @examples 
#' # Calculate relative adaptiveness values (ws) for E. coli K-12
#' eco.ws <- get.ws(tRNA=ecolik12$trna, sking=1)
#' 
#' # Calculate tAI for a set of 49 E. coli K-12 coding genes
#' eco.tai <- get.tai(ecolik12$m[,-33], eco.ws)
#' 
#' # Plot tAI vs. effective number of codons (Nc)
#' plot(eco.tai, ecolik12$w$Nc, xlab="tAI", ylab="Nc")
#' 
#' @export
get.tai <- function(x,w) {

  w = log(w)              #calculate log of w
  n = apply(x,1,'*',w)    #multiply each row of by the weights
  n = t(n)                #transpose
  n = apply(n,1,sum)      #sum rows
  L = apply(x,1,sum)      #get each ORF length
  tAI = exp(n/L)          #get tai
  return(tAI)
}

############################################################################
# After calculating tAI, we can then calculate S, i.e. the correlation
# between tAI and g(GC3s) - Nc (aka nc.adj)
############################################################################
#' Correlation between tAI and Nc adjusted
#' @author Mario dos Reis
#' @export
get.s <- function(tAI, nc, gc3) {

  nc.ad <- nc.adj(nc, gc3)
  ts <- cor(tAI, nc.ad, use='p')
  return(ts)
}

#############################################################################
# Statistical test for tAI
#############################################################################
#' Monte Carlo test of correlation between tAI and Nc adjusted
#' @author Mario dos Reis
#' @export
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
#' Nc vs. GC3s
#' 
#' Calculates the expected Nc value of a gene for a given GC content at the
#' third codon positions.
#' 
#' @param x a vector of GC contents at third codon positions
#' 
#' @details The relationship between Nc and GC content at third positions, x, 
#' is given by \deqn{Nc = -6 + x + 34/(x^2 + (1.025 - x)^2).} This equation
#' follows dos Reis et al. (2004, see also Wright 1990 for the original).
#' 
#' @return A vector of Nc values for the given GC contents.
#' 
#' @author Mario dos Reis
#' 
#' @references 
#' Wright F. (1990) The 'effective number of codons' used in a gene. \emph{Gene}, 
#' \bold{87:} 23--9.
#' 
#' dos Reis M., Savva R., and Wernisch L. (2004) Solving the riddle of codon 
#' usage preferences: a test for translational selection. \emph{Nucleic Acids Res.},
#' \bold{32:} 5036--44.
#' 
#' @examples
#' curve(nc.f(x), xlab="GC3s content", ylab="Nc")
#' points(ecolik12$w$GC3s, ecolik12$w$Nc, pch=19)
#' 
#' @export
nc.f <- function(x) {
  -6 + x + 34/(x^2 + (1.025 - x)^2)
}

# plot.nc <- function(nc, gc3s) {
#   curve(nc.f, xlim = c(0, 1), ylim = c(20, 62),
#         xlab = "GC3s", ylab = "Nc")
#   points(nc ~ gc3s, pch = '+', cex = 0.5)
# }
