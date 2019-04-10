# Bulk Segregant Analysis.

Our python program will help you find loci under selection with a bulk segregant approach and produce publication-quality figures using a simple command in Terminal. Although the programs were written and tested on sequences from the two-spotted spider mite, Tetranychus urticae, the python program is designed to work for related projects for which data sets and input files are available. We emphasize user-friendliness so do no hesitate to email Andre (a.kurlovs@gmail.com) if you have any questions. Bulk segregant analysis and this program is described in the following publication: 
- Kurlovs, A. H., R., Snoeck, Kosterlitz, O., Van Leeuwen, T., and Clark, R. M. Some flashy title. Under review (for a preprint, see bioRxiv ########; doi: some_magic_url)

---

## Installation
- Linux-based command line (Terminal). If you have MacOS, you just need to open Terminal. If you have Windows 10, you can run Linux command line as well – [check out this useful article](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
- Make sure [python](https://www.python.org/downloads/) 2 or 3 is installed.
- Install [numpy](https://scipy.org/install.html) 
- Install [Matplotlib](https://matplotlib.org/users/installing.html)

---

## Basic Run

### Sample experimental design
Suppose you are interested in finding the genomic location of insecticide (I) resistance loci. You individually cross ten insects from a resistant population (I_R) with a sensitive population (I_S), and expand the resulting populations. After several generations, you subject a subset of the population to the pesticide. The surviving mites are the target, resistant population (Selected). For each replicate, you then individually extract and sequence DNA from each bulk population; selected populations (Selected), unselected populations (Control), and both parents (I_R and I_S).
After you sequence each sample and acquire the fastq files, map each sample to the reference genome (recommend BWA; [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/)), and predict variants (recommend using GATK; [Genome Analysis Tool Kit](https://software.broadinstitute.org/gatk/)) to get a variant call format (VCF) file.

### Usage
Put the program in the working directory with your generated. To run the basic command, use the following template:  
```
/Users/Say_My_Name/My_BSA/bsa_code.py  -v /Users/Say_My_Name/Variant_Calls/my_stuff.vcf -psel I_R -pcon I_S -osel Selected_Offspring_1,Selected_Offspring_2,Selected_Offspring_3,Selected_Offspring_4, Selected_Offspring_5 -ocon Control_Offspring_1,Control_Offspring_2,Control_Offspring_3, Control_Offspring_4, Control_Offspring_5 -o /Users/Say_My_Name/My_BSA/Outfiles
```
To check all the options, run:
```
/Users/Say_My_Name/My_BSA/bsa_code.py –h
```
You may notice that you can provide single parental strains, and at the same time have multiple replicates. However, if you used multiple parental strains, they should also be separated by a comma and come in the same order as their respective offspring. Order matters!

This basic run will ourput the I_R allele frequency difference between each pair of selected and control offspring groups in sliding windows of 75kb (default option) that move across the genome in increments of 5kb (default option). It will output a number of files and plots. Note that the window size is based on the two-spotted spider mite (*Tetranychus urticae*) genome (~90Mb in size) and will need to be changed if your genome size is different. For corn (Zea mays), we used a window size of 5Mb and slide of 500kb (e.g.,`–w 5000000 –s 500000)`.

The basic run will output the following directories (within the major output directory you define as –o):

`/BSA_output` will have files from the sliding window analysis and statistics (see below).
`/BSA_plots` will by default have three plots:

`/BSA_plots/BSA_average_plot.pdf` plots a single line that represents the average of all selected-control offspring pairings.

`/BSA_plots/BSA_comb_plot.pdf` plot each selected control pairing individually.
`/BSA_plots/BSA_sep_plot.pdf` plots each selected and control sample separately in two different colors – one for selected offspring groups and one for the controls.

### Statistics

Our package performs a simulation (which is a type of permutation that shifts data around without changing the order) that determines whether your BSA peaks likely arose due to selection or genetic drift (please see [Wybouw et. al. 2019](https://doi.org/10.1534/genetics.118.301803) for explanation of this method). We use a false discovery rate (FDR) of 0.05 (default option), but it can also be defined by the user with the –sig flag (e.g., `–sig 0.01` for FDR of 0.01). The simulation relies on the power of replication – the more replicates you have, the better. It works by combining your selected samples with the respective control samples. If you want to run simulations on your data, add the –perm flag (e.g., `perm –i 10000`) to your command above. 

---

## Defaults

Maybe add a section here to talk about the assumptions of the data structure and design?

---

## Special Cases

### Unpaired Data

It often happens that (as used in our example) selected and control offspring are not paired and thus do not represent true replicates (please see Wybouw et. al. 2019 for example). Use `–u` for unpaired data. The order you put in is the order that will be used to pair them for plotting in BSA_average_plot.pdf. The crucial distinction is how permutations are performed. Because data are unpaired, every pairing has to be tested. Therefore, you may want to allow more time to process the permutations in this event, especially if you have many samples. By default, every sample combination will be permuted, which is the factorial for the number of samples. This means 120 in the example above, meaning that `–u –perm 10000` will perform a total of 1,200,000 permutations. 
With more permutations, this process will get very computationally intensive. Multiprocessing is incorporated into the code, and by default, the program will use all available processing cores on your machine. You can also specify how many cores you want use with the `–n` flag. To further reduce processing time, you can select how many random combinations you want to permute by using the –comb flag. For example, `–u –perm 10000 –n 30 –comb 60` will randomly select 60 selected-control offspring group combinations, and permute each 10k times using 30 processing cores. If you have many replicates (e.g., 10), the number of combinations rises quickly (e.g., 3.5 million). It is highly advised that you choose the number of desired combinations, e.g., `–comb 120`, if you have a large number of replicates.

### One parental strain

If you only have one parental strain sequenced, use the –pmaj flag. The alleles coming from the unknown parental strain will be inferred. For instance:
```
/Users/Say_My_Name/My_BSA/bsa_code.py  -v /Users/Say_My_Name/Variant_Calls/my_stuff.vcf -pmaj I_R -osel Selected_Offspring_1,Selected_Offspring_2,Selected_Offspring_3,Selected_Offspring_4, Selected_Offspring_5 -ocon Control_Offspring_1,Control_Offspring_2,Control_Offspring_3, Control_Offspring_4, Control_Offspring_5 -o /Users/Say_My_Name/My_BSA/Outfiles
```
Neither parental strain
If you have no parental information, do not specify the parents, and the program will perform BSA based on allele frequency change. Note that this method takes absolute value and is unidirectional. Important: if you have replicates, this method only works with paired data. Or rather, it will treat unpaired data as is does paired data.

### Heterozygous parents

The program’s default setting is to assume that both parental strains are homozygous. This will only focus on loci that are fixed for different alleles in both parents. However, if you have substantial heterozygosity in your parental strains, it may be worth trying our genotype inference methods. The two methods available are `–het` for when both parents are diploid and heterozygous and `–hpd` in which case you use a haploid male (from a potentially outbred/heterozygous strain).

–het 

This method only requires one of the parents to be fixed at a given locus, and it infers selected parent allele frequency in the experimental/control offspring groups by using presumed frequency in the F1 offspring. For instance, say your Parent I-R has A/G at a certain locus while parent I-S has A/A. If your offspring has both A and G alleles, I-R frequency will be calculated as freq(G)+0.5*freq(A). WARNING: Depending on your experimental design, this inference method may not be applicable at all. The design must be a _____

-hpd Haplodiploid_Male_Parent

This method has been created to improve inference for haplodiploid species. If a haploid male involved in a cross comes from a heterozygous strain, the allele that gets passed on to offspring can be inferred for each experimental/control offspring group. Using the example described above, if the I-R parent is a male from a heterozygous strain, and the offspring have both A and G alleles (with the A allele at less than 0.95 frequency by default, which can be changed using `–mac`, which is an acronym for “major allele cutoff”; the same was applied in the `–het` option), it is assumed that both alleles got passed on to offspring. This method also has shortcomings. The most glaring one is that an allele reaching frequency of over 0.95 can be in a region that is being swept to fixation because of selection. However, it is worth noting that the default method, which assumes homozygosity in both parental strains, would filter this allele out because it is a segregating site in one parent (A/G). 

---

## Plotting Customization

### Masking

While the actual BSA peak should have relatively smooth rise and fall, regions of misassembly in the genome can produce sudden sharp peaks. If it is reasonable to assume such regions are in fact misassembled, we provide an option to mask the region. The masking file should be created in a text editor (TextWrangler, Visual Studio Code, Notepad++, etc. Do not use programs like Word as they add special characters). The file should have tab-separated columns with chromosome (or scaffold name), beginning, and end position (bp) of the masking region. For instance,
```
chromosome_1	145000	155000
chromosome_4	1230000	1256000
```
After making the file, place in the current working directory for this project, e.g., 
`-mask /Users/Say_My_Name/My_BSA/masking_file.txt`

### Tickmarks and spacing

When the plot is automatically generated, the program attempts to space the tick marks on the x-axis based on the genome size, while the y-ticks are automatically spaced 0.1 apart. You can change these options. This can be handy for zooming in on a specific section of the genome (below). You can customize your plot with the following options (defaults shown below):
```
-xstep 0
-ystep 0.1
-xticks 17
-yticks 0
-xminor 5
-yminor 2
```
The “step” flags control how far apart the ticks are placed. The –xstep flag is set to 0 because the default plotting process will instead look a reasonable way to space the ticks (while aiming to make the total number of specified ticks, here it is set to 17 by `–xticks 17`) on the x-axis. Since the x-axis is longer, it will also include minor ticks that appear five times as frequently on the plot as the major ticks. You can control all of these options. However, one important thing to remember is that `–xstep` and `–ystep` will override `–xticks` and `–yticks`. In other words, the program won’t use `–xticks` or `–yticks` if you specify –xstep or `–ystep`. This means that since –ystep is provided by default, you should specify `–ystep 0` if you want your `–yticks` option to be used. It is important to point out that the automatic plotting methods should not be expected to always produce the best possible results. It may be helpful to play around with the options to generate what you think would be the best plot.

### Zoom

Sometimes it may be handy to zoom in on particular parts of the genome. To do so, you can provide a zoom file where each tab-delimited line (see MASKING above) specified the region you want to zoom in on. Add `-z /Users/Say_My_Name/My_BSA/zoom_file.txt` to your command. Note that since the default tick placing is 0.1, you should specify another tick spacing pattern if you want this to work well. For example, `–ystep 0 –yticks 10`. We also recommend that you use round numbers in your zoom file to avoid long decimals show up in the labels.

### Colors

The default color scheme is “jet”. For matplotlib color schemes, refer to <https://matplotlib.org/examples/color/colormaps_reference.html>.

The program takes schemes by default, but if you specify –col custom,
you can choose your own colors. Matplotlib has quite a few colors available. The github page has a pdf file with all the available colors (which I adapted from this post on Stack Overflow
<https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib>.)

If you want to choose your own color scheme, use e.g.,
`-col custom,red,blue,royalblue,orange,green,darkviolet,indianred`

`/BSA_plots/BSA_sep_plot.pdf` will have red for selected and blue for the unselected groups
`/BSA_plots/BSA_comb_plot.pdf` will use royalblue, orange, green, darkviolet, and indianred for your five replicates.

Alternatively, you can use hexadecimal color designations if you specify the `-hex` flag. Hexadecimals are to be used without the hash marks, e.g.,

`-col custom,B80F0A,111E6C,C29600,FC9C00,FAF00F`

### Independent plotting

The package allows users to plot BSA files of their choice and the statistical cutoffs associated with them. This can come in handy when you are analyzing three separate experiments with several replicates, yet want to produce a plot that summarizes all of your findings in the same place. This is accomplished by using `–plot` and `–permplot`. You can list your BSA sliding window files using the former and permutation files using the latter. For instance:
```
-plot /Users/Say_My_Name/My_BSA/Outfiles_experiment1/BSA_output/selected_average.txt, /Users/Say_My_Name/My_BSA/Outfiles_experiment2/BSA_output/selected_average.txt,
/Users/Say_My_Name/My_BSA/Outfiles_experiment3/BSA_output/selected_average.txt

-permplot
/Users/Say_My_Name/My_BSA/Outfiles_experiment1/BSA_output/permutations.txt, /Users/Say_My_Name/My_BSA/Outfiles_experiment2/BSA_output/permutations.txt,
/Users/Say_My_Name/My_BSA/Outfiles_experiment3/BSA_output/permutations.txt

-o /Users/Say_My_Name/My_BSA/New_outfiles
```
You do not need to input anything else to get the plot. However, and this is crucial – you do need to copy/paste an info_files directory into `/Users/Say_My_Name/My_BSA/New_outfiles` from one of your BSA runs. The program will automatically look for `/Users/Say_My_Name/My_BSA/New_outfiles/info_files/chrom_file.txt` to shade the chromosomes on your plot.

The colors will work the same in the sense that `–col` Accent will select that color scheme, while `–col custom,green,blue,violet` will use the three colors for your samples.


