Pedigree Simulator
==================
Program to simulate pedigrees containing arbitrary numbers of generations with
specified numbers of branches in each generation. The method can use
sex-specific genetic maps and randomly assigns the sex of each parent when
using such maps. Currently, all pedigrees must be symmetric and begin with the
oldest generation of founders that produce offspring that are either: full
siblings, half-siblings, or (in the third generation) double cousins. All
subsequent generations produce only full sibling descendents (although we aim to
make this more flexible).

Basic usage:

    ./ped-sim -d <in.dat> -m <map file> -i <in.vcf> -o <out_prefix>

The simulator produces four output files: `[out_prefix].vcf`, `[out_prefix].bp`,
`[out_prefix].fam`, and `[out_prefix].log`. Descriptions of each of these input
and output files appear below. Run `ped-sim` without arguments to see a full
listing of options; the non-required options are described at the end of this
document.

<!--TODO: give example; want to make small VCF with HapMap samples for this-->

------------------------------------------------------

Compiling
---------

To compile the ped-sim program, most Linux/Unix-based users should simply be
able to type

    make

Other setups may require editing of the Makefile or alternate means of
compiling.

<!--TODO: ideally distribute the binary-->

------------------------------------------------------

Dat file
--------

The dat file describes the pedigree structure(s) to be simulated. An example
dat file is given in `example/full_half_1st_2nd_cousin.dat` and a full
description of this file appears below.

The first line of a pedigree specification contains four columns:

    type #copies #generations name

Here, `type` is either `full`, `half`, or `double`, corresponding to simulating (more details given below):

* `full`: the samples in generation 2 are to be full siblings of each other
* `half`: the samples in generation 2 are to be half-siblings of each other
* `double`: the samples in generation 3 are to be double cousins of each other

`#copies` gives the number of replicated simulations of the given pedigree
structure to produce. These pedigrees will all have the same structure but will
descend from different founders, with different randomized sex assignments of
the founders, and so be independent.

`#generations` indicate the total number of generations in the pedigree.

`name` gives the name of the pedigree, which must be unique for each pedigree
structure in a given simulation run (i.e., a given dat file). The simulator
uses this to generate the simulated individuals' sample ids (details in the
"Sample ids for simulated individuals" section below).

After this first line, the dat file lists simulation details corresponding to
various generations in the pedigree. Each such line has the following format:

    generation# #samples_to_print [#branches]

`generation#` gives an integer value for the pedigree generation number. This
value can range from 1 (the earliest generation) to the total number of
generations included in the pedigree, as listed on the first line of the
specification (`#generations` just above).

`#samples_to_print` indicates how many samples the simulator should print from
each branch in the indicated generation. Note that the number of samples to
print does not change the number of branches: it merely allows for more
offspring to be simulated in the branches. Members of each branch have the same
parents, so multiple members of a given branch will be full siblings of each
other. As described further next, the branch count controls how many people in
that generation reproduce to form the next generation.

`#branches` is an optional field. In all generations but the last, exactly one
member of each branch reproduces with a single founder to produce at least
one or more offspring in a branch in the next generation. By default, for
`full` and `half` type pedigrees, generation 2 contains two branches, each with
one offspring of the parents in generation 1. (The members of these two
branches are either full siblings of each other or are half-siblings
corresponding to `full` and `half` type pedigrees, respectively.) And by
default, every generation includes the same number of branches as the previous
generation. The branch number in every generation must be greater than or equal
to the immediately previous generation. If the branch number is greater than
the previous generation, it must be an integer multiple of the previous value.
When the branch number increases by a multiple of *n*, branch *i* generates
branches *n\*(i-1)+1* through *n\*i*, where *i* ranges from 1 to \#branches
(from the previous generation).

<!--TODO: talk about #branches in double-->

The following dat entry simulates five three-generation pedigrees, with two
branches (the default) in generation 2, both branches containing four children
for a total of eight full siblings in the second generation. As there are two
branches in the second generation, there are only two individuals in generation
2 that reproduce. These individuals each have two children in the two branches
within generation 3, yielding a total of four grandchildren.

    full 5 3
    2 4
    3 2

If data for the third generation samples is all that is needed, the following
is equivalent to the above, since `#samples_to_print` defaults to 0:

    full 5 3
    3 2

The following dat entry produces four branches in generation 2, each of which
produce children (grandchildren), but does not print any samples from
generation 2.

    full 5 3
    2 0 4
    3 2

<!--TODO: probably give an error when the branch count in the last generation
is not equal to the previous generation-->

Thus this gives a total of 4\*2 = 8 grandchildren with four sets of branches
that each contain two full siblings. The four sets are first cousins of one
another.

The dat file can contain any number of distinct pedigree structures, which the
program simulates separately using unique founders for all pedigrees (including
multiple copies of the same structure).

Comments are allowed on a line by themselves beginning with #.

### Example dat file: `example/full_half_1st_2nd_cousins.dat`

The repository includes a dat file (in `example/full_half_1st_2nd_cousins.dat`)
that simulates two different but closely related pedigree structures. The first
is a `full` type pedigree so that the children in generation 2 are full siblings
of each other. The file does not list entries for either generation 1 or 2, and
therefore produces the default of two branches in generation 2 and does not
print the simulated individuals from either of these generations. The
specification calls for four branches in generation 3, so both children in
generation 2 produce two children in generation 3, each belong to its own
branch. Generation 4 is the last generation and produces and prints 1 child per
existing branch, for a total of four individuals produced. Because the four
branches in generation 3 included two sets of full siblings, two pairs of the
four samples in generation 4 are first cousins. The other pairs are second
cousins, with their most recent common ancestors in generation 1.

The second pedigree structure differs from the first only in that it is a `half`
type pedigree. This results in the two children that belong to the two
(default) branches in generation 2 being half-siblings rather than full
siblings. As a consequence, in generation 4, there are two sets of half-first
cousins and the remaining individuals are half-second cousins.

------------------------------------------------------

Map file
--------

The genetic map file contains three columns for a sex-averaged map or four
columns if using male and female maps. The format is:

    chromosome physical_position map_position0 [map_position1]

The chromosomes are expected to be listed in the same order as they appear in
the input VCF file, with the physical positions in increasing order. These
chromosomes names must match the names in the input VCF file.

`map_position0` is genetic position in centiMorgans, and should either be the
sex-averaged genetic position if using only one map, or should be the male
genetic position if using two maps. When using only one map, the simulator
samples all recombinations from that one map and does not distinguish (nor
assign) male and female parents.

`map_position1` is likewise a genetic position in centiMorgans and should
correspond to the female genetic position if given.

A high resolution sex-specific genetic map is available [here](https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination),
and is described in [Bhérer, et al. 2017](http://dx.doi.org/10.1038/ncomms14994).
To generate a map file in the format the simulator requires with both male and
female genetic positions, the following works (run through the bash shell):

    printf "#chr\tpos\tmale_cM\tfemale_cM\n" > refined_mf.simmap
    for chr in {1..22}; do
      paste male_chr$chr.txt female_chr$chr.txt | awk -v OFS="\t" 'NR > 1 && $2 == $6 {print $1,$2,$4,$8}' | sed 's/^chr//' >> refined_mf.simmap;
    done

This generates a file called `refined_mf.simmap` that can be passed to the
simulator.

------------------------------------------------------

Input VCF file
--------------

All founders in the simulated pedigrees are a randomly sampled input individual
from the input VCF file. This VCF must contain phased data for all individuals,
with no missing data for any site. As most phasers automatically impute missing
data, the latter requirement should be easily met.

------------------------------------------------------

Output VCF file
---------------

The output VCF contains the simulated individuals, including only those samples
requested to be printed according to the dat file. For any generation in which
there is a request to print one or more samples, the simulator prints any
founders in that generation as well as the non-founders. See below for a
description of the sample ids for the simulated individuals.

------------------------------------------------------

Output BP file
--------------

The break points (BP) file lists information about the composite haplotypes of
each sample. All founders have a unique numerical id for each of their two
haplotypes, starting from 0 and ranging to 2\**F*-1, where *F* is the number of
founders across all simulated pedigrees. Within the BP file, there are two
lines for every sample that is to be printed (according to the dat file). Each
line begins with the sample id of the simulated individual (described below),
the sex of that person, either s0 for male or s1 for female, a designation of
the haplotype that that line describes, either h0 or h1, and then a variable
number of segments for each simulated chromosome.

<!--TODO: describe remainder of BP file-->

------------------------------------------------------

Output fam file
---------------

The simulator produces a PLINK format fam file with the pedigree structures
simulated. This fam file contains all generated samples, including those that
are not requested to be printed in the dat file. This enables the relationships
between all samples to be determined from the fam file alone.

------------------------------------------------------

Output log file
---------------

Information about the simulation run appears in the log file and is a copy of
what is printed to the console during execution. Notably this includes the
random seed used for a given simulation. Supplying the same input files with
the same random seed (assignable with the `--seed` option) will produce the
same simulation results.

------------------------------------------------------

Sample ids for simulated individuals
------------------------------------

The simulated individuals' sample ids have the format
`[name][#]_g[#]-b[#]-i[#]`. Here, `[name]` is the pedigree name given in the
dat file. The first number `[#]` is the copy number of the pedigree which
ranges from 1 to the number of copies requested in the dat file. The `g[#]`
portion of the id gives the generation number of the individual, which ranges
from 1 to the total number of generations in the pedigree. `b[#]` gives the
branch number the sample is contained in in the indicated generation; this
ranges from 1 to the total number of branches in that generation. Finally,
`i[#]` gives the individual number in the given branch and generation. This
ranges from 0 to the total number of samples requested to be simulated in the
generation (inclusive). Thus, for each branch, the output contains the number
of samples to be printed plus (in most cases) individual 0. For all generations
but the last, individual 1 and a founder individual 0 are the parents of the
branches they produce in the next generation. (The last generation does not
reproduce so does not include any founders or samples labeled individual 0.)
The only exception to the rule that individual 0 is a founder is in generation 2
in double cousin simulations: indivuals 0 and 1 in both branches are siblings of
each other and marry an individual from the other branch (who descend from
distinct founding parents).

------------------------------------------------------

Extra notes
-----------

When simulating with sex-specific maps, it is necessary to include data for all
chromosomes in one run. This is because sex is assigned randomly, but only once
per run. Thus, to maintain consistency of the sex of each individual in a
given pedigree, all chromosomes need to be simulated in the same run.

------------------------------------------------------

Other optional arguments
------------------------

### Specifying random seed: `--seed <#>`

The `--seed <#>` option enables specification of the random seed to be used.
Without this option, the simulator generates a random seed using the current
time (including microseconds).

### Genotyping error rate: `--err_rate <#>`

To more accurately mimic real data, the simulator introduces genotyping errors
at a specified rate, defaulting to 1e-3. Set this value to 0 to keep the
allelic values identical to those in the founder haplotypes (from the input
data).

**Note: only pedigree samples have genotyping errors introduced;
`--retain\_extra` samples maintain their original calls**

### Rate of opposite homozygote errors: `--err_hom_rate <#>`

SNP array genotype calling works by clustering allele intensities among a set
of samples. So if an individual is truly homozygous, its intensities are more
likely to fall in either the correct cluster or the heterozygous cluster, with
a lower probability of being called homozygous for the opposite allele.  While
we are unaware of a study that looks at error rates by "true" genotype class in
SNP array data, the `--err_hom_rate` option provides the ability to produce
different rates of errors for genotypes that are truly homozygous. The default
rate for generating an erroneous genotype that is homozygous for the opposite
alleles relative to the truth is 0.1. That is, when a homozygous genotype is
set to an erronous value, 10% of the time it is set homozygous for the opposite
allele, and 90% of the time is heterozygous. For equal rates of both these
classes, set the rate for this option to .5.

### Missingness rate: `--miss_rate <#>`

As real data includes missingness, the simulator introduces missing genotype
calls at a rate specified by this parameter, with a default of 5e-3. Set this
value to 0 for no missing genotypes.

**Note: only pedigree samples have sites set to missing; `--retain\_extra`
samples maintain their original calls**

### Maintaining phase in output: `--keep_phase`

By default the simualtor produces a VCF that does not contain phase information.
The `--keep_phase` option will instead generate a VCF that mainitains the
phase of all samples.

### Retaining extra input samples: `--retain_extra <#>`

The simulator uses samples from the input VCF as founder individuals, and will
exit if more founders are needed than available in the input VCF. If requested
using `--retain_extra <#>` the program will also print a specified number of
input samples that were not used as founders in the simulations. If the number
is less than 0 (e.g., `--retain_extra -1`), the simulator prints all unused
input samples. If the value is greater than 0, say 100, but fewer than this
number of unused samples exist, the simulator prints all the available samples.
When the requested number to print is less than the number available, the
simulator randomly selects the samples to print from among all that were not
used as founders.
