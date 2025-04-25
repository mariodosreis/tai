# The tRNA adaptation index

tAI is an R package (and two perl scripts) for the analysis of codon usage in
DNA coding sequences. It implements the tRNA adaptation index (tAI) described in
dos Reis et al. (2003, 2004). The tAI package is distributed WITHOUT WARRANTY
under the terms of the GNU General Public License. See the file LICENSE for
details.

This package is basically the same I used in my papers on tAI, with just a few
minor modifications to make it a little friendlier. If you have any doubts,
questions, bug reports, etc., please contact me at:

Mario dos Reis

mariodosreis@gmail.com

School of Biological and Behavioural Sciences  
Queen Mary University of London


INSTALLATION
----------------------------------------------------------------------

tAI is now available from CRAN. From the R prompt:

~~~R
install.packages("tAI")
~~~

Alternatively, (assumming you have the devtools package installed) you can also install tAI from GitHub:

~~~R
devtools::install_github("mariodosreis/tai")
~~~

And that's it! The package has been installed!

QUICK TUTORIAL
----------------------------------------------------------------------

You are assumed to know how to use a shell or command prompt, as well as some
basic knowledge of R. `$` denotes the shell prompt (so don't type it in!).

The tRNA adaptation index (tAI) measures the degree of co-adaptation between a
coding sequence and the tRNA pool of an organism (dos Reis 2003, 2004). File
`inst/extdata/ecolik12.ffn` contains 49 coding sequences extracted from the
genome of *Escherichia coli* K-12. We can use the files in the package to
calculate tAI for each one of the genes contained in this file. First we need to
calculate the frequencies of the 61 coding codons for every sequence. Copy the
files inside the `inst/extdata` directory into a new suitable directory in your
system. Go into your new directory and type (note that perl must be installed in
your system and in your path):

~~~
$ perl codonM ecolik12.ffn ecolik12.m
~~~

The file `ecolik12.m` contains the output of the codonM script: a matrix of
codon frequencies per ORF. It should look like:

~~~
11      19      10      13      11      10      6       9       ...
6       4       3       6       0       6       1       3       ...
13      11      5       12      0       2       3       5       ...
0       0       2       0       1       1       1       0       ...
8       8       2       6       0       4       2       1       ...
...
~~~

Each row represents one ORF or gene (in our case, there should be 49 rows in the
file `ecolik12.m`), and the columns represent each one of the 61 coding codons,
arranged in this fashion:

~~~
1      TTT      14     CTT      30     ATT      46     GTT
2      TTC      15     CTC      31     ATC      47     GTC
3      TTA      16     CTA      32     ATA      48     GTA
4      TTG      17     CTG      33     ATG      49     GTG

5      TCT      18     CCT      34     ACT      50     GCT
6      TCC      19     CCC      35     ACC      51     GCC
7      TCA      20     CCA      36     ACA      52     GCA
8      TCG      21     CCG      37     ACG      53     GCG

9      TAT      22     CAT      38     AAT      54     GAT
10     TAC      23     CAC      39     AAC      55     GAC
-      -        24     CAA      40     AAA      56     GAA
-      -        25     CAG      41     AAG      57     GAG

11     TGT      26     CGT      42     AGT      58     GGT
12     TGC      27     CGC      43     AGC      59     GGC
-      -        28     CGA      44     AGA      60     GGA
13     TGG      29     CGG      45     AGG      61     GGG
~~~

Notice that STOP codons have been excluded. Also, codonM ignores the first codon
in every sequence, this is because it is always a Methionine codon (even if its
not coded by the canonical ATG). The codons above follow the TCAG ordering. The
standard genetic code ordered this way is

~~~
AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
~~~

Now start R in the same directory you have been working. Type:

~~~R
require("tAI")
~~~

This will load the necessary functions into R. The file `ecolik12.trna` contains
the gene copy number of every kind of tRNA in the *E. coli* K-12 genome. We need
this information in order to calculate tAI:

~~~R
eco.trna <- scan("ecolik12.trna")
~~~

This file contains 64 rows, corresponding to the anticodon complements of each
tRNA species (e.g. if the tRNA anticodon is GAA, the complement is TTC). The
tRNAs  are ordered according to their anticodon complement, in the same order as
in codonM's output as indicated above, but with added STOP codons. STOP codons
are at positions 11, 12 and 15 as in the standard genetic code diagram above.
Now we can calculate the relative adaptiveness values for each codon in the E.
coli genes:

~~~R
eco.ws <- get.ws(tRNA=eco.trna, sking=1)
~~~

Now, lets load the output of codonM into R:

~~~R
eco.m <- matrix(scan("ecolik12.m"), ncol=61, byrow=TRUE)
~~~

We will ignore Methionine codons in our analysis:

~~~R
eco.m <- eco.m[,-33]
~~~

Now we can finally calculate tAI:

~~~R
eco.tai <- get.tai(eco.m, eco.ws)
hist(eco.tai)
~~~

The last command will plot an histogram of the tAI values. Highly expressed
genes present high tAI values (> 0.4), which means that their codon usage
resembles the genomic structure of tRNA genes.

Now a good question arises, how much of the total codon usage of these genes (in
terms of Nc) is due to adaptation to the tRNA gene pool? In order to answer this
question, we need to calculate Nc for every sequence in `ecolik12.ffn`. You
should have codonW (or other codon usage package able to calculate Nc)
installed. If you do have codonW, and it is in your path, you can use codonZ to
quickly and efficiently compute a set of codon usage statistics. At the shell
type:

~~~
$ perl codonZ ecolik12.ffn ecolik12.w
~~~

If you don't have codonW installed, file `ecolik12.w` is already provided in the
`inst/extdata` directory. Now from the R prompt:

~~~R
df <- read.table("ecolik12.w", header=TRUE)
~~~

If you are using the codonW version distributed by me, the above command should
work fine, if you are working with the original distribution, you should type
instead:

~~~R
df <- read.table("ecolik12.w", header=TRUE, na.strings = "*****")
~~~

Lets plot the relationship between tAI and Nc:

~~~R
plot(eco.tai ~ df$Nc)
~~~

As you can see, genes with very low Nc values (highly biased), correlate with
high tAI values (highly co-adapted to the tRNA gene pool).

~~~R
cor(eco.tai, df$Nc, use="p")
~~~

You should get a value of –0.9100338

Formally, we can calculate the correlation between tAI, and the corrected Nc,
f(GC3s) – Nc, we call this correlation S, because it reflects the intensity of
translational selection acting on our sample of genes:

~~~R
eco.s <- get.s(eco.tai, df$Nc, df$GC3s)
~~~

You should get a value of 0.9065442

For more details, read the references!

REFERENCES
----------------------------------------------------------------------

1. dos Reis M, Savva R, and Wernisch L. (2003) Unexpected correlations between
gene expression and codon usage bias from microarray data for the whole
Escherichia coli K-12 genome. *Nucleic Acids Research*, **31:** 6976–6985.
[DOI: 10.1093/nar/gkg897](https://dx.doi.org/10.1093/nar/gkg897)

2. dos Reis M, Wernisch L, and Savva R. (2004) Solving the riddle of codon usage
preferences: a test for translational selection. *Nucleic Acids Research*,
**32:** 5036–5044.
[DOI: 10.1093/nar/gkh834](https://dx.doi.org/10.1093/nar/gkh834)
