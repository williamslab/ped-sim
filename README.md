Pedigree Simulator
==================
Program to simulate pedigrees containing arbitrary numbers of generations with
specified numbers of branches in each generation. The method can use
sex-specific genetic maps and randomly assigns the sex of each parent when
using such maps. Currently, all pedigrees must be symmetric and begin with the
oldest generation of founders that produce offspring that are either: full
siblings, half-siblings, or (in the third generation) double cousins.

Usage:

    ./ped-sim [in.dat] [map file] [in.vcf] [out_base] <random seed>

The simulator produces four output files: `[out_base].vcf`, `[out_base].bp`,
`[out_base].fam`, and `[out_base].log`. Descriptions of each of these files
appear below.

<!--TODO: give example-->

------------------------------------------------------

Compiling
---------

To compile the ped-sim program, most Linux/Unix-based users should simply be
able to type

    make

Other setups may require editing of the Makefile or alternate means of
compiling.

------------------------------------------------------

Dat file
--------

The dat file describes the pedigree structure(s) to be simulated. An example
dat file is given in `example/full_half_1st_2nd_cousin.dat` and this file is
described more fully later.

The first line of a pedigree specification contains three columns:

    type #copies #generations name

Here, `type` is either `full`, `half`, or `double`, corresponding to simulating (more details given below):

* `full`: the samples in generation 2 are to be full siblings of each other
* `half`: the samples in generation 2 are to be half-siblings of each other
* `double`: the samples in generation 3 are to be double cousins of each other

`#copies` gives the number of replicated simulations of the given pedigree
structure to produce. These pedigrees will all have the same structure but will
descend from different founders and so be independent.

`#generations` indicate the total number of generations in the pedigree.

`name` gives the name of the pedigree, which must be unique for each pedigree
structure in a given simulation run (i.e., a given dat file). This is used for
generating the simulated individuals' sample ids (see below).

After this first line, the dat file lists simulation details corresponding to
various generations in the pedigree. Each such line has the following format:

    generation# #samples_to_print <#branches>

`generation#` gives an integer value for the pedigree generation number. This
value can range from 1 (the earliest generation) to the total number of
generations included in the pedigree, as listed on the first line of the
specification (`#generations` just above).

`#samples_to_print` indicates how many samples should be printed from each
branch in the indicated generation. Note that the number of samples to print
does not change the number of branches: it merely allows for more offspring to
be simulated in the extant branches. Each branch has the same parents, so
multiple members of a given branch will be full siblings of each other. As
described further next, the branch count controls how many people in that
generation reproduce to form the next generation.

`#branches` is an optional field. In all generations but the last, exactly one
member of each branch reproduces with a single non-founder to produce at least
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
branches *n\*i* through *n\*(i+1)-1*, where *i* ranges from 0 to \#branches-1.

<!--TODO: talk about #branches in double-->

The following simulates five three-generation pedigrees, with two branches
(the default) in generation 2, both branches containing four children for
a total of eight full siblings in the second generation. As there are two
branches in the second generation, there are only two individuals in generation
2 that reproduce. These children each have two children in each branch,
yielding a total of four grandchildren.

    full 5 3
    2 4
    3 2

If data for the third generation samples is all that is needed to be printed,
the following is equivalent to the above, since `#samples_to_print` defaults to
0:

    full 5 3
    3 2

The following dat entry produces four branches in generation 2, each of which
has their own children (grandchildren) but does not print any samples from
generation 2.

    full 5 3
    2 0 4
    3 2

<!--TODO: probably give an error when the branch count in the last generation
is not equal to the previous generation-->

Thus this gives a total of 4\*2 = 8 grandchildren with four sets of branches
that each contain two full siblings. The four sets are first cousins of one
another.

The dat file can contain any number of distinct pedigree structures, each of
the program simulates separately, using unique founders for all pedigrees
(including multiple copies of the same structure).

Comments are allowed on a line by themselves beginning with #.

### Example dat file: `example/full_half_1st_2nd_cousins.dat`

The repository includes a dat file (in `example/full_half_1st_2nd_cousins.dat`)
that simulates two different but closely related pedigree structures. The first
is a `full` type pedigree so that the children in generation 2 are full siblings
of each other. The file does not list entries for either generation 1 or 2, and
therefore produces the default of two branches in generation 2. The
specification calls for four branches in generation 3, so both children in
generation 2 produce two children in generation 3, each belong to its own
branch. Generation 4 is the last generation and produces and prints 1 child per
existing branch, for a total of four individuals produced. Because the four
branches in generation 3 included two sets of full siblings, two pairs of these
four samples are first cousins. The other pairs are second cousins, with their
most recent common ancestors in generation 1.

The second pedigree structure differs from the first only in that it is a `half`
type pedigree. This results in the two children that belong to the two (default)
branches in generation 2 are half-siblings of each other rather than full
siblings. In consequence, in generation 4, there are two sets of half-first
cousins and the remaining individuals are half-second cousins.

------------------------------------------------------

Map file
--------

The genetic map file contains three columns for a sex-averaged map or four
columns if using male and female maps. The format is:

    chromosome physical_position map_position0 <map_position1>

The chromosomes are expected to be listed in the same order as they appear in
the input VCF file, with the physical positions in increasing order.

`map_position0` is genetic position in centiMorgans, and should either be the
sex-averaged genetic position if using only one map, or should be the male
genetic position if using two maps. When using only one map, the simulator
samples all recombinations from that one map and does not distinguish (nor
assign) male and female parents.

`map_position1` is likewise a genetic position in centiMorgans and should
correspond to the female genetic position.

<!--TODO: link to 23andMe map-->

------------------------------------------------------

Input VCF file
--------------

All founders in the simulated pedigrees are a randomly sampled input individual
from the input VCF file. The VCF must contain phased data for all individuals.

------------------------------------------------------

Output VCF file
---------------

The output VCF contains the simulated individuals, including only those samples
to be printed according to the dat file. For any generation in which there is a
request to print one or more samples, the simulator prints any founders in that
generation as well as the non-founders. The file also includes all individuals
from the input VCF that were not assigned as a founder. Thus, if the input VCF
contains 100 samples and the simulator only requires 10 founders to perform the
simulation, the output VCF will contain all the simulated individuals as well
as the 90 samples that were not used as founders for the simulation. Founder
assignment is randomized, so different runs using the same input VCF and dat
file will produce differ

See below for a description of the sample ids for the simulated individuals.

------------------------------------------------------

Output BP file
--------------

The break points (BP) file lists information about which haplotypes each sample
carries. All founders have a unique numerical id for each of their two
haplotypes, starting from 0 and ranging to 2\**F*-1, where *F* is the number of
founders across all simulated pedigrees. Within the BP file, there are two
lines for every sample that is to be printed (according to the dat file). Each
line beings with the sample id of the simulated individual (as described
below), the sex of that person, either s0 for male or s1 for female, the
haplotype that that line describes, h0 or h1, and then a variable number of
segments for each simulated chromosome.

<!--TODO: describe remainder of BP file-->

------------------------------------------------------

Output fam file
---------------

The simulator produces a PLINK format fam file with the pedigree structure of
the individuals. This fam file includes all generated samples, including those
that are not requested to be printed in the dat file. This enables the
relationships between all samples to be determined from the fam file alone.

------------------------------------------------------

Output log file
---------------

Information about the simulation run, a copy of what is printed to the console
during execution, appears in the log file. Notably this includes the random seed
used for a given simulation. Supplying the same input files along with a given
random seed to the simulator will produce the same simulation data.

------------------------------------------------------

Sample ids for simulated individuals
------------------------------------

The simulated individuals' sample ids have the format
`[name][#]_g[#]_b[#]_i[#]`. Here, `[name]` is the pedigree name given in the
dat file. The first number `[#]` is the copy number of the pedigree which
ranges from 1 to the number of copies requested in the dat file. The `g[#]`
portion of the id gives the generation number of the individual, which ranges
from 1 to the total number of generations in the pedigree. `b[#]` gives the
branch number the sample is contained in in the indicated generation; this
ranges from 1 to the total number of branches in that generation. Finally,
`i[#]` gives the individual number in the given branch and generation. This
ranges from 0 to the total number of samples requested to be simulated in the
generation (inclusive). Thus there can be as many as the number of samples to
be printed plus an additional sample: sample 0 is generally a founder individual
(with the exception being for generation 2 in double cousin simulations).

------------------------------------------------------

Extra notes
-----------

When simulating with sex-specific maps, it is necessary to include data for all
chromosomes in one run. This is because sex is assigned randomly, but only once
per run. Thus, to maintain consistency of the sex of each individual in a
given pedigree, all chromosomes need to be simulated in the same run.
