# Haplo
Haplo is an implementation of the algorithm described here: http://dx.doi.org/10.5048/BIO-C.2016.4 
It is partial in that population substructure has not yet been included, but all other features are included.

Haplo performs a full backwards coalescent simulation of a chromosome or section of chromosome, for arbitrarily large populations and (almost) arbitrarily long timescales. It includes the possibility to simulate different numbers of males and females, different mating behaviours, different (or changing) mutation rates, different recombination rates. Full details of all the options are given below.

### COMPILE, SET-UP, AND RUN
---

**Compile** using `hpc/build.sh` to create an executable `Haplo`. 
(Make sure you have the latest version of g++ and make).
The script will tell you where to find the built executable.
Copy the executable to wherever you will run it.  

**Edit** input files, for example named `inputs.txt`. Examples can be found in `hpc/sample_inputs`.
Typically you will need to create a population-structure-history file, and add its name to the input file. 
You can also create a temporary (scratch) directory and/or an output directory and add them to the input file.

**Run** with `./Haplo inputs.txt`.
The code makes checkpoints periodically. If the program gets stopped by causes beyond its control (if it does not crash),
restart it with `./Haplo inputs.txt --restart`
Take care restarting with a changed input file. For example, if you change `temporary-file-prefix` input, it may not be able to find the saved data it needs. But restarting can sometimes be useful even if the run has completed. For example, one can rerun the final analysis with a different set of analysis parameters, without having to re-run the whole simulation. 

To stop a simulation cleanly before it has completed, create a file named `stop` (for example, type `touch stop`) in the working directory. If the stop file is present when the code next reaches a suitable stopping point, it will exit cleanly.


### INPUT FILE
---
The input file uses a version of TOML composed of key-value pairs in the form `key = value`.
The `#` character indicates that any following text on a line is a comment to be ignored.
Lines may be continued onto the next line by beginning the next line with whitespace.

##### GENETICS 
`chromosome_type`   - value can be `autosome` or `x` or `y` or `mito`
 `chromosome_length` - chromosome length in base pairs
`recombination_rate` - per generation per nucleotide
`mutation_rate`      - per generation per nucleotide
`mutation_rate_history_file` - if empty use constant mutation_rate, if file given see details below

##### PRIMORDIAL STRUCTURE (ID THEORY)
`primordial_diversity` - the heterozygosity : probability that any given nucleotide is heterozygous
`primordial_block_length` - length of primordial blocks in base pairs
`primordial_probability_dimorphic` - the proportion of blocks that have two variants
`primordial_probability_tetramorphic` - the proportion of blocks that have four variants

##### POPULATION / MATING
`population_structure_history_file` - file containing demographic history - see below for details
`fertility_parameter_alpha` - any number from 0 to infinity
    - infinity (or say 10^10) means mothers have random numbers of children
    - small numbers mean some mothers are more likely to have children
    - 0 means one mother has all children
`mating_parameter_beta`- any number from 0 to infinity
    - infinity (or say 10^10) means mothers mate randomly with fathers
    - small numbers mean mothers are more likely to have children with same father
    - 0 means each mother has children with only one father

##### SAMPLING AND CULLING
`random_seed` - change this to create a new sample history with the same parameters
`population_sample` - the number of individuals to be sampled in extant generation
`maximum_blocks` - the number of haplotype blocks that can be created (limit for compute efficiency)
`use_mutation_loci` - `true`: store positions of mutations, `false`: approximate them for efficiency
`cull_non_parents` - `true/false` remove non-parents from further consideration?
`cull_nonancestral_parents` - `true/false` remove parents with no ancestral material?
`hide_nonancestral_blocks` - `true/false` remove/hide blocks with no ancestral material?

##### IMPLEMENTATION 
`verbose_logging` - `true/false` extra logging comments
`use_memory_unloading` - `true/false` unload memory contents to temporary files?
`temporary_file_prefix` - prefix or path for naming or placing temporary files
`output_file_prefix` - prefix or path for naming or placing output files

##### ANALYSIS / OUTPUT
`analysis_datapoints_maximum` - maximum resolution for all output data distributions
`analysis_linkage_distance_maximum` - maximum distance to calculate linkage
`analysis_do_linkage_stats` - `true/false` option to do or not-do expensive linkage calculation
`analysis_linkage_minimum_frequency` - value from 0 to 0.5
    -exclude alleles with frequency smaller than this from linkage calculation
`output_final_SNP_sites` - `true/false` option to create vcf like output of all SNPs (or not)

### POPULATION-STRUCTURE-HISTORY FILE
---
This file contains a series of triples: a generation number, and then the number of males and the number of females. The size of each subpopulation is interpolated between each point. The generation numbers are counted back in time, but they are listed forwards in time from the founding generation to the present. Therefore the generation numbers must be decreasing. For example, for a simulation of 5000 generations with a constant population of 1000 men and 1500 women per generation:
```
5000 1000 1500
0    1000 1500
```
For a population that slowly grows from 2000 to 4000 and then plateaus at 4000:
``` 
5000 1000 1000
4000 2000 2000
0    2000 2000
```

### MUTATION-RATE-HISTORY FILE
---
If the filename is not set, the constant `mutation_rate` parameter (see above) is used instead. This file is very similar to the population-structure-history file except that it has a series of pairs: a generation number and the mutation rate at that generation. Mutation rates are again linearly interpolated between points. 
For example, for mutation rate that falls linearly from 2 x 10^-8 to 1 x 10^-8 over 2000 generations:
```
2000 2e-8
0 1e-8
```
---
---
