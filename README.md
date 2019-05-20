Our program will help you find loci under selection with a bulked segregant approach and produce publication-quality figures. We emphasize user-friendliness so do not hesitate to email Andre (a.kurlovs@gmail.com) if you have any questions or concerns. This program is described in the following publication:
- Kurlovs, A. H., Snoeck, S., Kosterlitz, O., Van Leeuwen, T., and Clark, R. M. Some flashy title. Under review (for a preprint, see bioRxiv ########; doi: some_magic_url)

---

# Requirements
- Linux-based command line (Terminal). If you have MacOS, you just need to open Terminal. If you have Windows 10, you can run Linux command line as well – [check out this useful article](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
- Make sure [Python](https://www.python.org/downloads/) 3 is installed
- Install Python packages [numpy](https://scipy.org/install.html) and [Matplotlib](https://matplotlib.org/users/installing.html)

---

# Sample Experimental Design and Basic Run

## Input files and experimental design
Suppose you are interested in finding the genomic location of insecticide (I) resistance loci. First, you individually cross five insects from a resistant population (I_R) with a sensitive population (I_S), and expand the resulting segregating populations. After several generations, subject a subset of each of the five replicates to the pesticide. The surviving mites become your selected populations (referred to as Selected_Offspring in the example below) while the total populations are the unselected (referred to as Control_Offspring in the example below). You then extract and sequence DNA from each selected and unselected population as well as from both parental strains (I_R and I_S). More information and suggestions on experimental design for BSA studies are available in the publication specified above.
After you sequence each sample and get the fastq files containing the sequence reads, you map each sample to the reference genome (recommend BWA; [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/)), and predict variants (recommend using GATK; [Genome Analysis Tool Kit](https://software.broadinstitute.org/gatk/)) to get a variant call format (VCF) file, which will be used as the input file for program.

## Usage
You get your VCF file and you are ready to find your pesticide resistance locus. To run the basic command, use the following template:  
```
python /Users/Say_My_Name/My_BSA/bsa_code.py \
-v /Users/Say_My_Name/Variant_Calls/my_stuff.vcf \
-psel I_R \
-pcon I_S \
-osel Selected_Offspring_1,Selected_Offspring_2,Selected_Offspring_3,Selected_Offspring_4,Selected_Offspring_5 \
-ocon Control_Offspring_1,Control_Offspring_2,Control_Offspring_3,Control_Offspring_4,Control_Offspring_5 \
-o /Users/Say_My_Name/My_BSA/Outfiles
```
To read more about these options and check all of them, run:
```
/Users/Say_My_Name/My_BSA/bsa_code.py –h
```
The basic command is best suited for a situation in which both of your parental strains have been inbred and sequenced. This method only considers homozygous loci, which might limit your power. Please see the [Two inbred parental strains not present](#Two-inbred-parental-strains-not-present) section to see options that deal with other scenarios.

You may notice that you can provide single parental strains and at the same time have multiple replicates. However, if you use multiple parental strains, they should also be separated by a comma and come in the same order as their respective offspring. Order matters! 

This basic run will determine the I_R allele frequency difference between each pair of selected and control offspring groups in sliding windows of 75kb (default option) that move across the genome in increments of 5kb (default option). It will output a number of files and plots. Note that the window size is based on the two-spotted spider mite (*Tetranychus urticae*) genome (~90Mb in size) and will need to be changed if your genome size is different. For corn (Zea mays), we use a window size of 5Mb and a slide of 500kb (e.g., `–w 5000000 –s 500000`).

The basic run will output the following directories (within the major output directory you define as `–o`):

`/BSA_output` will have files from the sliding window analysis and statistics (see below).
`/BSA_plots` will by default have three plots:

`/BSA_plots/BSA_average_plot.pdf` plots a single line that represents the average of all selected-control offspring pairings.

`/BSA_plots/BSA_comb_plot.pdf` plots each selected control pairing individually.
`/BSA_plots/BSA_sep_plot.pdf` plots each selected and control sample separately in two different colors – one for the selected offspring groups and one for the controls.

---

# Statistics - Is it likely selection or drift? 

## Overview
Our package performs a simulation (which is a type of permutation that shifts allele frequency data around without changing the order) that determines whether your BSA peaks likely arose due to selection or genetic drift (please see [Wybouw, Kosterlitz, et. al. 2019](https://doi.org/10.1534/genetics.118.301803) for an explanation of this method). We use a false discovery rate (FDR) of 0.05 (default option), but it can also be defined by the user with the –sig flag (e.g., `–sig 0.01` for FDR of 0.01). The simulation relies on the power of replication – the more replicates you have, the better. It works by combining your selected samples with the respective control samples. If you used different parental strains in your BSA, you will need to run the code individually for each set of parents for this to work correctly. If you want to run simulations on your data, add the –perm flag (e.g., `perm –i 10000`) to your command above. This method assumes that your selected and control methods are paired. Refer to the section below ([Treating unpaired data](#treating-unpaired-data)) if otherwise. 

## Treating unpaired data
It often happens that selected and control offspring are not paired and thus do not represent true replicates (please see the tomato selection experiment in [Wybouw, Kosterlitz, et. al. 2019,](https://doi.org/10.1534/genetics.118.301803) for an example). Use `–u` for unpaired data. The order you put in is the order that will be used to pair them for plotting in `BSA_average_plot.pdf`. The crucial distinction is how permutations are performed. Because data are unpaired, every potential pairing has to be tested. Therefore, you may want to allow more time to process the permutations in this event, especially if you have many samples. By default, every sample combination will be permuted, which is the factorial of the number of samples. This means that if you have five replicates, `–u –perm 10000` will perform a total of 1,200,000 permutations. 
With more permutations, this process will get very computationally intensive. Multiprocessing is incorporated into the code, and by default, the program will use all available processing cores on your machine. You can also specify how many cores you want to use with the `–n` flag. To further reduce processing time, you can select how many random combinations you want to permute by using the `–comb` flag. For example, `–u –perm 10000 –n 30 –comb 60` will randomly select 60 selected-control offspring group combinations, and permute each 10k times using 30 processing cores. If you have many replicates (e.g., 10), the number of combinations rises quickly (with 10 replicates, to over 3.5 million). It is thus highly advised that you choose the number of desired combinations, e.g., `–comb 120`, if you have a large number of replicates.

---

# Two inbred parental strains not present

We recommend using inbred parental strains and sequencing them in the course of your experiment. However, we realize that this is not always possible. The methods described in this section might remedy your situation. 


## Haplodiploid male parent
As this method was developed for the two-spotted spider mite, we offer an option in which individual haploid male parents from a heterozygous strain were crossed to a female from an inbred strain. This option assumes that your samples are paired as different males would be used in each cross. It does not matter whether a strain with the trait of interest or a strain without it the haplodiploid one. Use the `–hpd' flag for the haplodiplod male parental strain, e.g.,
```
/Users/Say_My_Name/My_BSA/bsa_code.py \
-v /Users/Say_My_Name/Variant_Calls/my_stuff.vcf \
-psel I_R \
-pcon I_S \
-hpd I_R \
-osel Selected_Offspring_1,Selected_Offspring_2,Selected_Offspring_3,Selected_Offspring_4,Selected_Offspring_5 \
-ocon Control_Offspring_1,Control_Offspring_2,Control_Offspring_3,Control_Offspring_4,Control_Offspring_5 \
-o /Users/Say_My_Name/My_BSA/Outfiles
```

## One inbred parental strain
If you only have one of the two parental strains inbred and sequenced, use the `–pmaj` flag. The alleles coming from the unknown parental strain will be inferred. Like the '-hpd' method described above, this one treats your samples as if they are paired. It is also built for the parent with the trait being available -- otherwise, your BSA peak will point downwards.
```
/Users/Say_My_Name/My_BSA/bsa_code.py \
-v /Users/Say_My_Name/Variant_Calls/my_stuff.vcf \
-pmaj I_R \
-osel Selected_Offspring_1,Selected_Offspring_2,Selected_Offspring_3,Selected_Offspring_4,Selected_Offspring_5 \
-ocon Control_Offspring_1,Control_Offspring_2,Control_Offspring_3,Control_Offspring_4,Control_Offspring_5 \
-o /Users/Say_My_Name/My_BSA/Outfiles
```

## No genetic information on parental strains
If you have no parental information, do not specify the parents, and the program will perform BSA based on the absolute allele frequency difference between the selected and the unselected offspring. This also treats data as paired, and since it has to combine the pairs for each allele, it makes permutations not possible.

---

# Additional information

## Troubleshooting

In a sliding window analysis, the minimum number of SNPs that have to be in a window to be considered is by default set as window_size*0.0005. Given that the default window size is 75kb, at last 38 SNPs have to be present in the window for it to be considered. If your parental strains have few SNPs compared to the reference genome, this setting might be too stringent. Use the `–m` flag followed the minimum number of SNPs of your choice to change this parameter. 

## Masking

While the actual BSA peak should have relatively smooth rise and fall, regions of misassembly in the genome can produce sudden sharp peaks. If it is reasonable to assume such regions are in fact misassembled, we provide an option to mask the region. The masking file should be created in a text editor (TextWrangler, Visual Studio Code, Notepad++, etc. Do not use programs like Word as they add special characters). The file should have tab-separated columns with chromosome (or scaffold name), beginning, and end position (bp) of the region to mask. For instance,
```
chromosome_1	145000	155000
chromosome_4	1230000	1256000
```
After making the file, you specify its location, e.g., 
`-mask /Users/Say_My_Name/My_BSA/masking_file.txt`

# Plotting parameters

## Tick marks and spacing
When the plot is automatically generated, the program attempts to space the tick marks on the x-axis based on the genome size, while the y-ticks are automatically spaced 0.1 apart. You can change these options. This can come particularly handy when you are zooming in on a specific section of the genome (below). You can customize your plot with the following options (defaults shown below):
```
-xstep 0
-ystep 0.1
-xticks 17
-yticks 0
-xminor 5
-yminor 2
```
The “step” flags control how far apart the ticks are placed. The `–xstep` flag is set to 0 because the default plotting process will instead look a reasonable way to space the ticks (while aiming to make the total number of specified ticks; here it is set to 17 by `–xticks 17`) on the x-axis. Since the x-axis is longer, it will also include minor ticks that appear five times as frequently on the plot as the major ticks. You can control all of these options. However, one important thing to remember is that `–xstep` and `–ystep` will override `–xticks` and `–yticks`. In other words, the program won’t use `–xticks` or `–yticks` if you specify `–xstep` or `–ystep`. This means that since `–ystep` is provided by default, you should specify `–ystep 0` if you want your `–yticks` option to be used. It is important to point out that the automatic plotting methods should not be expected to always produce the best possible results. It may be helpful to adjust options to generate what you think would be the best plot.

## Zoom
Sometimes it may be handy to zoom in on particular parts of the genome. To do so, you can provide a zoom file where each tab-delimited line (see [Masking](#masking) above) specified the region you want to zoom in on. Add `-z /Users/Say_My_Name/My_BSA/zoom_file.txt` to your command. Note that since the default tick placing is 0.1, you should specify another tick spacing pattern if you want this to work well. For example, `–ystep 0 –yticks 10`. We also recommend that you use round numbers in your zoom file to avoid long decimals show up in the labels.

## Colors
The default color scheme is “jet”. For Matplotlib color schemes, refer to <https://matplotlib.org/examples/color/colormaps_reference.html>.

The program takes schemes by default, but if you specify `–col` custom,
you can choose your own colors. The github page has a pdf file with all the available color names (which we adapted from this post on Stack Overflow
<https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib>). Alternatively, you can use hexadecimal color designations.

If you want to choose your own color scheme, use e.g.,
`-col custom,red,blue,royalblue,orange,#00FF00,darkviolet,indianred`

`/BSA_plots/BSA_sep_plot.pdf` will have red for selected and blue for the unselected groups
`/BSA_plots/BSA_comb_plot.pdf` will use royalblue, orange, #00FF00 (which is hex for lime green), darkviolet, and indianred for your five replicates.

# Independent plotting

The package allows users to plot BSA files of their choice and the statistical cutoffs associated with them. This can come in handy when you are analyzing three separate experiments with several replicates, yet want to produce a plot that summarizes all of your findings in the same place. This is accomplished by using `–plot` and `–permplot`. You can list your BSA sliding window files using the former and permutation files using the latter. For instance:
```
-plot /Users/Say_My_Name/My_BSA/Outfiles_experiment1/BSA_output/selected_average.txt,
/Users/Say_My_Name/My_BSA/Outfiles_experiment2/BSA_output/selected_average.txt,
/Users/Say_My_Name/My_BSA/Outfiles_experiment3/BSA_output/selected_average.txt

-permplot
/Users/Say_My_Name/My_BSA/Outfiles_experiment1/BSA_output/permutations.txt,
/Users/Say_My_Name/My_BSA/Outfiles_experiment2/BSA_output/permutations.txt,
/Users/Say_My_Name/My_BSA/Outfiles_experiment3/BSA_output/permutations.txt

-o /Users/Say_My_Name/My_BSA/New_outfiles
```
You do not need to input anything else to get the plot. However, and this is crucial – you do need to copy/paste an info_files directory into `/Users/Say_My_Name/My_BSA/New_outfiles` from one of your BSA runs. The program will automatically look for `/Users/Say_My_Name/My_BSA/New_outfiles/info_files/chrom_file.txt` to shade the chromosomes on your plot.

The colors will work the same in the sense that `–col` Accent will select that color scheme, while `–col custom,green,blue,violet` will use the three colors for your samples.



