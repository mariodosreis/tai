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
#' 
#' The adjusted Nc is f(gc3s) - Nc
#' 
#' @param nc a vector of length n with the effective number of codons for genes
#' @param gc3 a vector of length n with corresponding GC composition at third codon positions
#' 
#' @details The adjusted Nc is calculated as described in dos Reis et al. (2004).
#' 
#' @references
#' dos Reis M., Savva R., and Wernisch L. (2004) Solving the riddle of codon 
#' usage preferences: a test for translational selection. \emph{Nucleic Acids Res.},
#' \bold{32:} 5036--44.
#' 
#' @seealso \code{\link{nc.f}} for the function used to calculate f(gc3s)
#' 
#' @author Mario dos Reis
#' 
#' @examples 
#' 
#' eco.ncadj <- nc.adj(ecolik12$w$Nc, ecolik12$w$GC3s)
#' plot(eco.ncadj ~ ecolik12$w$Nc, xlab="Nc", ylab="Nc adjusted")
#' 
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
#' Calculates the relative adaptiveness values of codons based on the number 
#' of tRNA genes.
#' 
#' @param tRNA a vector of length 64 with tRNA gene copy numbers
#' @param s a vector of length 9 with selection penalties for codons
#' @param sking a vector of length 1 indicating the superkingdom
#' 
#' @details The relative adaptiveness values are calculated as described in 
#' dos Reis et al. (2003, 2004). If \code{s = NULL}, the s values are set to 
#' the optimised values of dos Reis et al. (2004). \code{sking} indicates the 
#' superkingdom, with 0 indicating Eukaryota, and 1 Prokaryota.
#' 
#' @return A vector of length 60 of relative adaptiveness values.
#' 
#' @author Mario dos Reis
#' 
#' @examples 
#' eco.ws <- get.ws(tRNA=ecolik12$trna, sking=1)
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
#' @export
get.ws <- function(tRNA,      # tRNA gene copy number
                   s = NULL,  # selective constraints
                   sking)     # super kingdom: 0-eukaryota, 1-prokaryota
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
#' @param x an n by 60 matrix of codon frequencies for n open reading frames.
#' @param w a vector of length 60 of relative adaptiveness values for codons.
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
#' @seealso \code{\link{get.ws}}
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
#' 
#' Calculates the correlation between tAI and Nc (adjusted for GC content at
#' third codon positions).
#' 
#' @param tAI a vector of length n with tAI values for genes
#' @param nc a vector of length n with Nc values for genes
#' @param gc3 a vector of length n with GC content at third codon positions for genes
#' 
#' @author Mario dos Reis
#' @importFrom stats cor
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
#' 
#' Calculates the p-value (using a Monte Carlo or randomisation test) that the 
#' correlation (the S value) between tAI and the adjusted Nc for a set of 
#' genes is different from zero.
#' 
#' @param m a k by 60 matrix of codon frequencies for k genes
#' @param ws vector of length 60 of relative adaptiveness values of codons
#' @param nc vector of length k of Nc values for genes
#' @param gc3s vector of length k of GC content at third codon position for genes
#' @param ts.obs vector of length 1 with observed correlation between tAI and Nc adjusted for the k genes
#' @param samp.size a vector of length 1 with the number of genes to be sampled from m (see details)
#' @param n the number of permutations of ws in the randomisation test
#' 
#' @details The Monte Carlo test is described in dos Reis et al. (2004). When
#' working with complete genomes, matrix \code{m} can have a very large number 
#' of rows (large k). In this case it may be advisable to choose \code{samp.size}
#' < k to speed up the computation.
#' 
#' @return A list with elements \code{p.value}, the p-value for the test, and 
#' \code{ts.simulated}, a vector of length \code{n} with the simulated 
#' correlations between tAI and adjusted Nc.
#' 
#' @examples 
#' eco.ws <- get.ws(tRNA=ecolik12$trna, sking=1)
#' eco.tai <- get.tai(ecolik12$m[,-33], eco.ws)
#' ts.obs <- get.s(eco.tai, ecolik12$w$Nc, ecolik12$w$GC3s)
#' 
#' # The S-value (dos Reis et al. 2004):
#' ts.obs # [1] 0.9065442
#' 
#' # There seems to be a high correlation between tAI and Nc adjusted for
#' # the 49 genes in ecolik12$m. Is the correlation statistically significant?
#' ts.mc <- ts.test(ecolik12$m[,-33], eco.ws, ecolik12$w$Nc, ecolik12$w$GC3s, 
#'                  ts.obs, samp.size=dim(ecolik12$m)[1])
#' # The p-value is zero:
#' ts.mc$p.value # [1] 0
#' 
#' # Histogram of simulated S-values:
#' hist(ts.mc$ts.simulated, n=50, xlab = "Simulated S values", 
#'      xlim=c(min(ts.mc$ts.simulated), ts.obs))
#' # Add the observed S-value as a red vertical line:
#' abline(v=ts.obs, col="red")
#' 
#' @references 
#' dos Reis M., Savva R., and Wernisch L. (2004) Solving the riddle of codon 
#' usage preferences: a test for translational selection. \emph{Nucleic Acids Res.},
#' \bold{32:} 5036--44.
#' 
#' @author Mario dos Reis
#' @export
ts.test <- function(m, ws, nc, gc3s, ts.obs, samp.size, n=1000) {

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

  # p-values is:
  p.value <- sum(ts > ts.obs)

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
#' @details Without selection on codon bias, the expected value of Nc as a 
#' function of GC content at third positions, x, is given by 
#' \deqn{f(x) = -6 + x + 34/(x^2 + (1.025 - x)^2).} This equation
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

#############################################################################
# Data documentation
#############################################################################
#' E. coli K-12 codon bias and tRNA numbers
#' 
#' A list with elements \code{trna}, a vector of length 64 of tRNA gene copy numbers 
#' in the Escherichia coli K-12 genome, \code{w}, a data frame with some codon bias 
#' statistics for 49 E. coli K-12 coding genes, and \code{m}, a 49 by 61 matrix of 
#' codon frequencies for the 49 genes in question.
#'
#' @author Mario dos Reis
#' 
#' @examples 
#' # 87 tRNA genes in the E. coli K-12 genome:
#' sum(ecolik12$trna)
#' 
#' # Two copies are isoacceptors for Phe, with anticodon GAA (codon TTC)
#' ecolik12$trna[2]
#' 
#' # ecolik12$w, a data frame with codon bias statistics
#' names(ecolik12$w)
#' 
#' # Effective number of codons vs. gene length (in codons)
#' plot(ecolik12$w$Nc, ecolik12$w$L_aa, xlab="Nc", ylab="Gene length")
"ecolik12"