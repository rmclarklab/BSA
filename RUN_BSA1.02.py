#!/usr/bin/env python

"""This program performs BSA and associated statistics"""

import sys
import re
import argparse
import math
import multiprocessing
import os
import random
import subprocess

from decimal import Decimal
from itertools import permutations
from numpy import percentile

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
from matplotlib import pyplot as plt

N_CPU = multiprocessing.cpu_count()
rcParams['font.sans-serif'] = 'Arial'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

PARSER = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

# vcf file and out-directory
PARSER.add_argument("-v", "--vcf", required=False,
                    help="VCF file with BSA parents and offspring")
PARSER.add_argument("-o", "--outdir", required=True,
                    help="Directory where output files will be written")
# parent and offspring linermation and individual variant filters
PARSER.add_argument("-psel", "--selected_parent", required=False,
                    help="Parent with trait of interest "
                    "if several crosses, insert them in order, commas,no_spaces")
PARSER.add_argument("-pcon", "--control_parent", required=False,
                    help="Parent without trait of interest "
                    "if several crosses, insert them in order, commas,no_spaces")
PARSER.add_argument("-pmaj", "--major_parent", required=False,
                    help="Known parent when the other's genotype is not known "
                    "if several crosses, insert them in order, commas,no_spaces")
PARSER.add_argument("-hpd", "--haplodiploid", required=False, default=None,
                    help="Haplodiploid male parent")
PARSER.add_argument("-osel", "--selected_offspring", required=False,
                    help="Offspring with trait of interest "
                    "if several crosses, insert them in order, commas,no_spaces")
PARSER.add_argument("-ocon", "--control_offspring", required=False,
                    help="Offspring without trait of interest "
                    "if several crosses, insert them in order, commas,no_spaces")
PARSER.add_argument("-mac", "--mac", required=False, default=0.95,
                    help="Major allele cutoff; "
                    "locus will not be considered segregating if higher")
PARSER.add_argument("-over", "--coverage_over", required=False, default=1.50,
                    help="Maximum coverage depth as a multiple of genome-wide average")
PARSER.add_argument("-under", "--coverage_under", required=False, default=0.25,
                    help="Minimum coverage depth as a multiple of genome-wide average")
PARSER.add_argument("-qds", "--qds", required=False, default=2,
                    help="Mininum variant score")
PARSER.add_argument("-sor", "--sor", required=False, default=3,
                    help="Maximum strand bias score")
PARSER.add_argument("-mq", "--mps", required=False, default=50,
                    help="Minimum mean mapping quality score")
PARSER.add_argument("-mqrs", "--mqrs", required=False, default=-8,
                    help="Minimum mean quality rank sum score")
PARSER.add_argument("-rprs", "--rprs", required=False, default=-8,
                    help="Minimum read pos rank sum")
# General options for parsing the VCF file and for sliding window analysis
PARSER.add_argument("-b", "--binsize", required=False, default=100000,
                    help="Genomic size of bins for storing VCF linermation")
PARSER.add_argument("-w", "--window", required=False, default=75000,
                    help="Genomic length of sliding windows")
PARSER.add_argument("-s", "--slide", required=False, default=5000,
                    help="Genomic length by which windows slide across the genome")
PARSER.add_argument("-m", "--min_allele", required=False,
                    help="Minimum number of variants in a window")
PARSER.add_argument("-f", "--min_scaffold", required=False, default=500000,
                    help="Minimum chromsome/scaffold length")
#stuff for plotting
PARSER.add_argument("-xstep", "--xstep", required=False, default=0,
                    help="Distance between x-marks")
PARSER.add_argument("-ystep", "--ystep", required=False, default=0.1,
                    help="Distance between y-marks")
PARSER.add_argument("-xticks", "--xticks", required=False, default=17,
                    help="Number of ticks on x-axis")
PARSER.add_argument("-yticks", "--yticks", required=False, default=0,
                    help="Number of ticks on y-axis")
PARSER.add_argument("-xminor", "--xminor", required=False, default=5,
                    help="Number of minor x-ticks related to spacing of major ticks")
PARSER.add_argument("-yminor", "--yminor", required=False, default=2,
                    help="Number of minor y-ticks related to spacing of major ticks")
PARSER.add_argument("-col", "--colors", required=False,
                    default='jet',
                    help="Colors for plot")
PARSER.add_argument("-plot", "--plot", required=False, default=None,
                    help="Combined sel-control files separated by comma")
PARSER.add_argument("-permplot", "--permplot", required=False, default=None,
                    help="Permutation files separated by comma")
# this is for permutations
PARSER.add_argument("-perm", "--perm", required=False, default=0,
                    help="The number of permutations to perform")
PARSER.add_argument("-sig", "--significance", required=False, default=0.05,
                    help="Significance cutoff for permutation test")
PARSER.add_argument("-sigcol", "--sigcolor", required=False, default="red",
                    help="Color of line denoting significance on plot")
PARSER.add_argument("-u", "--unpaired", required=False, action="store_true",
                    help="Use if experimental and control groups unpaired; "
                    "by default the are paired")
PARSER.add_argument("-n", "--n_threads", required=False, default=N_CPU,
                    help="Number of threads; "
                         "only used when data are unpaired; "
                         "defaults to the number of processing core")
PARSER.add_argument("-comb", "--combinations", required=False, default=1,
                    help="Number of exp-control combinations "
                         "that will be permuted "
                         "the default stays 1 if paired"
                         "and factorial(ngroups) if not"
                         "if unpaired and any number above 1"
                         "it will use that number")
# additional options
PARSER.add_argument("-mask", "--masking_file", required=False, default=None,
                    help="File with genomic regions to mask"
                         "Chrom\tbeg\tend\n for masking")
PARSER.add_argument("-zoom", "--zoom_file", required=False, default=None,
                    help="File with regions where you want a zoomed in plot"
                         "Chrom\tbeg\tend\n for zoomed in plotting")
PARSER.add_argument("-vb", "--verbose", required=False, action="store_true",
                    help="Prints additional files")

########################## VARIABLES ############################

ARGDICT = {}

ARGIES = PARSER.parse_args()
ARGDICT["vcf"] = ARGIES.vcf
ARGDICT["outdir"] = ARGIES.outdir
ARGDICT["outdir1"] = ARGDICT["outdir"]+"/info_files"
ARGDICT["outdir2"] = ARGDICT["outdir"]+"/BSA_output"
ARGDICT["outdir3"] = ARGDICT["outdir"]+"/BSA_plots"

ALL_POPS = set()

if ARGIES.selected_parent:
    ARGDICT["selected_parent"] = ARGIES.selected_parent.split(",")
    ALL_POPS = ALL_POPS | set(ARGDICT["selected_parent"])
if ARGIES.control_parent:
    ARGDICT["control_parent"] = ARGIES.control_parent.split(",")
    ALL_POPS = ALL_POPS | set(ARGDICT["control_parent"])
if ARGIES.haplodiploid:
    ARGDICT["haplodiploid"] = ARGIES.haplodiploid
if ARGIES.major_parent:
    ARGDICT["major_parent"] = ARGIES.major_parent.split(",")
    ALL_POPS = ALL_POPS | set(ARGDICT["major_parent"])
if ARGIES.selected_offspring:
    ARGDICT["selected_offspring"] = ARGIES.selected_offspring.split(",")
    ALL_POPS = ALL_POPS | set(ARGDICT["selected_offspring"])
if ARGIES.control_offspring:
    ARGDICT["control_offspring"] = ARGIES.control_offspring.split(",")
    ALL_POPS = ALL_POPS | set(ARGDICT["control_offspring"])

ARGDICT["mac"] = float(ARGIES.mac)
ARGDICT["coverage_over"] = float(ARGIES.coverage_over)
ARGDICT["coverage_under"] = float(ARGIES.coverage_under)

ARGDICT["binsize"] = int(ARGIES.binsize)
ARGDICT["window"] = int(ARGIES.window)
ARGDICT["slide"] = int(ARGIES.slide)
if ARGIES.min_allele:
    ARGDICT["min_allele"] = int(ARGIES.min_allele)
else:
    ARGDICT["min_allele"] = ARGDICT["window"]*0.00050
ARGDICT["min_scaffold"] = int(ARGIES.min_scaffold)
ARGDICT["qds"] = float(ARGIES.qds)
ARGDICT["mps"] = float(ARGIES.mps)
ARGDICT["sor"] = float(ARGIES.sor)
ARGDICT["mqrs"] = float(ARGIES.mqrs)
ARGDICT["rprs"] = float(ARGIES.rprs)

ARGDICT["xstep"] = float(ARGIES.xstep)
ARGDICT["ystep"] = float(ARGIES.ystep)
ARGDICT["xticks"] = int(ARGIES.xticks)
ARGDICT["yticks"] = int(ARGIES.yticks)
ARGDICT["xminor"] = float(ARGIES.xminor)
ARGDICT["yminor"] = float(ARGIES.yminor)
if ARGIES.plot:
    ARGDICT["plot"] = ARGIES.plot.split(",")
if ARGIES.permplot:
    ARGDICT["permplot"] = ARGIES.permplot.split(",")
       
ARGDICT["color"] = ARGIES.colors.split(",")

if len(ARGDICT["color"]) == 1:
    ARGDICT["color"] = ARGDICT["color"][0]
else:
    ALL_COLORS_NAME = set([])
    for collie in matplotlib.colors.cnames.items():
        ALL_COLORS_NAME.add(str(collie[0]))
    INVALID_COLORS = []
    for col in ARGDICT["color"][1:]:
        if col[0] == "#":
            match_color = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', col)
            if not match_color:
                INVALID_COLORS.append(col)
        else:
            if col not in ALL_COLORS_NAME:
                INVALID_COLORS.append(col)
    if INVALID_COLORS:
        INVALID_STRING = ",".join(INVALID_COLORS)
        print("THE FOLLOWING COLORS ARE INVALID: %s"%(INVALID_STRING))
        sys.exit()
    if "selected_offspring" in ARGDICT:
        if len(ARGDICT["color"])-3 != len(ARGDICT["selected_offspring"]):
            print("INVALID NUMBER OF COLORS! NEED TO PROVIDE COLORS FOR: "
                "SELECTED, UNSELECTED, AND ALL THE REPLICATES")
            sys.exit()
    if "plot" in ARGDICT:
        if len(ARGDICT["color"])-1 != len(ARGDICT["plot"]):
            print("INVALID NUMBER OF COLORS!")
            sys.exit()

ARGDICT["perm"] = int(ARGIES.perm)
ARGDICT["sig"] = float(ARGIES.significance)
ARGDICT["sigcolor"] = ARGIES.sigcolor
if ARGIES.unpaired:
    ARGDICT["unpaired"] = ARGIES.unpaired
    ARGDICT["n_threads"] = int(ARGIES.n_threads)
    if int(ARGIES.combinations) > 1:
        ARGDICT["combinations"] = int(ARGIES.combinations)
    else:
        ARGDICT["combinations"] = math.factorial(len(ARGDICT["control_offspring"]))
else:
    ARGDICT["n_threads"] = 1
    ARGDICT["combinations"] = 1

if ARGIES.masking_file:
    ARGDICT["masking_file"] = ARGIES.masking_file
if ARGIES.verbose:
    ARGDICT["verbose"] = ARGIES.verbose
if ARGIES.zoom_file:
    ARGDICT["zoom_file"] = ARGIES.zoom_file

PREFIX = {0:'bp',
          1:'kb', # kilo
          2:'Mb', # mega
          3:'Gb', # giga
          4:'Tb', # tera
          5:'Pb', # peta
          6:'Eb', # exa
          7:'Zb', # zetta
          8:'Yb'  # yotta
         }

########################## FUNCTIONS ############################

def mean(vector):
    """Calculates mean"""
    list_vector = list(vector)
    sum_vector = sum(list_vector)
    len_vector = float(len(list_vector))
    mean_final = sum_vector/len_vector
    return(mean_final)

def error(message):
    """Prints error messages"""
    sys.exit(message)

def arange(start, end, step):
    """A range function that works for floats"""
    split_start = str(start).split(".")
    if len(split_start) > 1:
        dec_start = len(split_start[1])
    else:
        dec_start = 0
    split_step = str(step).split(".")
    if len(split_step) > 1:
        dec_step = len(split_step[1])
    else:
        dec_step = 0
    if dec_step > dec_start:
        decider = str(step)
    else:
        decider = str(start)
    final_list = []
    current = float(start)
    while current < end:
        current = int(current) if float(current) == int(current) else current
        final_list.append(current)
        current = current + step
        current = float(Decimal(current).quantize(Decimal(decider)))
    return(final_list)

def afill(start, end, ntries):
    """A function that fill evenly spaced values between two numbers"""
    step = (end-start)/float(ntries+1) if ntries > 0 else 0
    final_list = [float(start) + (i+1)*step for i in range(ntries)]
    return(final_list)

# this function processes scaffold linermation
def scale_dict(line_file):
    """Creates a dict with cumulative positions for each chromosome"""
    line_read = open(line_file, "r")
    line_dict = {}
    final_lenny = 0
    lenny = 0
    scaffolds = []
    for line in line_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        scaffolds.append(scaffy)
        final_lenny = final_lenny + lenny
        lenny = int(linetab[1])
        line_dict[scaffy] = final_lenny
    line_read.close()
    return(line_dict, scaffolds)

# this function processes scaffold linermation
def length_dict(line_file):
    """Creates a dict with chromosome lengths"""
    line_read = open(line_file, "r")
    line_dict = {}
    for line in line_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        lenny = int(linetab[1])
        line_dict[scaffy] = lenny
    line_read.close()
    return(line_dict)

def scale_ends(line_file):
    """Creates a dict with ends for each chromosome"""
    try:
        line_read = open(line_file, "r")
    except IOError:
        error("FILE NOT FOUND, LIKELY BECAUSE YOU ARE RUNNING "
              "THE PLOTTING CODE. MOVE THE info_files DIRECTORY "
              "INTO YOUR OUTDIR FOR THIS TO WORK")
    line_dict = {}
    final_lenny = 0
    lenny = 0
    scaffolds = []
    for line in line_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        scaffolds.append(scaffy)
        lenny = int(linetab[1])
        final_lenny = final_lenny + lenny
        line_dict[scaffy] = final_lenny
    line_read.close()
    return(line_dict, scaffolds)

def masker():
    """Specifies regions of genome to mask"""
    maskdict = {}
    if "masking_file" in ARGDICT:
        maskfile = open(ARGDICT["masking_file"])
        for line in maskfile:
            line = line.split("\t")
            try:
                int(line[1]) and int(line[2])
            except ValueError:
                error("CANNOT CONVERT MASKING COORDS TO INTEREGRS. "
                      "EXITING PROGRAM.")
            if line[0] not in maskdict:
                maskdict[line[0]] = set(range(int(line[1]), int(line[2])+1))
            maskdict[line[0]] = maskdict[line[0]] | set(range(int(line[1]), int(line[2])+1))
    scaffs = scale_dict(ARGDICT["outdir1"] + "/chrom_file.txt")[1]
    for scaff in scaffs:
        if scaff not in maskdict:
            maskdict[scaff] = set([-1])
    return(maskdict)

def verbosy(indict, contigs, outkeys):
    """Outputs all the variant/allele info"""
    for outkey in outkeys:
        outfile = open(ARGDICT["outdir1"] + "/%s.txt"%(outkey), "w")
        outfile.close()
    for contig in contigs:
        binnies = indict[contig]
        for binny in sorted(binnies):
            positions = indict[contig][binny]
            for position in sorted(positions):
                strains = indict[contig][binny][position]
                for strain in strains:
                    value = indict[contig][binny][position][strain]
                    outfile = open(ARGDICT["outdir1"] + "/%s.txt"%(strain), "a")
                    # should remove binny once i test this
                    outfile.write("%s\t%s\t%s\t%s\n"%(contig, position, binny, value))
                    outfile.close()
# COVERAGE
def coverage():
    """Calculates average coverage per strain/individual"""
    print("ITERATING OVER VCF FILE TO GET COVERAGE INFORMATION")
    cov = {}
    if not os.path.isdir("%s"%ARGDICT["outdir1"]):
        subprocess.call("mkdir %s"%(ARGDICT["outdir1"]), shell=True)
    outcov = open(ARGDICT["outdir1"]+"/coverageinfo.txt", "w")
    line_file_handle = open(ARGDICT["outdir1"]+"/chrom_file.txt", "w")
    the_right_stuff = set()
    openvcf = open(ARGDICT["vcf"])
    for line in openvcf:
        if line[0] == "#":
            if line.split("=")[0] == "##contig":
                line = line.split("=")
                vcf_scaff = line[2].split(",")[0]
                vcf_scaff_len = (line[3].split(">")[0])
                if int(vcf_scaff_len) > ARGDICT["min_scaffold"]:
                    the_right_stuff.add(vcf_scaff)
                    line_file_handle.write("%s\t%s\n"%(vcf_scaff, vcf_scaff_len))
            elif line[0:6] == "#CHROM":
                header_strains = line.strip().split("\t")[9:]
                not_in_vcf = ALL_POPS - set(header_strains)
                if len(not_in_vcf) > 0:
                    not_in_vcf_string = ",".join(not_in_vcf)
                    print("INVALID STRAIN NAMES: %s. EXITING PROGRAM"%(not_in_vcf_string))
                    sys.exit()
                right_strains = [st for st in header_strains if st in ALL_POPS]
                indexed_strains = [header_strains.index(st) for st in right_strains]
                for strain in right_strains:
                    cov[strain] = {}
                    cov[strain]["cov"] = 0.0
                    cov[strain]["total"] = 0.0
        else:
            line = line.replace("|", "/")
            line = line.strip().split("\t")
            info = line[9:]
            if (line[0] in the_right_stuff
                    and len(line[3]) == 1
                    and len(line[4].split(",")) == 1):
                for ix_strain in zip(indexed_strains, right_strains):
                    if "./" not in info[ix_strain[0]]:
                        reads = info[ix_strain[0]].split(":")[1].split(",")
                        parent_cov = sum([int(i) for i in reads])
                        cov[ix_strain[1]]["cov"] = cov[ix_strain[1]]["cov"]  + parent_cov
                        cov[ix_strain[1]]["total"] = cov[ix_strain[1]]["total"] + 1
    openvcf.close()
    line_file_handle.close()
    for strain in right_strains:
        avecov = cov[strain]["cov"]/cov[strain]["total"]
        outcov.write("%s\t%s\n"%(strain, avecov))
    outcov.close()

def spt_vcfcov(stringy):
    """Parses string with coverage info in VCF and gets total coverage"""
    return(sum([float(i) for i in
                stringy.split(":")[1].split(",")]))

def round_sig(xxx, sig):
    """Determines the number of sig figs"""
    if xxx:
        output = (round(xxx, sig-int(math.floor(
            math.log10(abs(xxx))))-1))
    else:
        output = 0
    return(output)

def inferred(indv):
    """Process offspring information when information only from a single parent"""
    exp_genosplito = set(indv[2].split(":")[0].split("/"))
    sample_scores = []
    if len(exp_genosplito) == 1:
        exp_allele = list(exp_genosplito)[0]
        #### now let's get the selected
        for ns in range(2):
            sample_ix = list(set([call for call in indv[ns].split(":")[0].split("/")
                                if call == exp_allele]))
            if sample_ix:
                sample_ix = int(sample_ix[0])
                sample_score = float(indv[ns].split(":")[1].split(",")
                                    [sample_ix])/spt_vcfcov(indv[ns])
            else:
                sample_score = 0.0
            sample_scores.append(sample_score)
        high_scores = [i for i in sample_scores if i >= ARGDICT["mac"]]
        if len(high_scores) < len(sample_scores):
            return(sample_scores)

def haplodiploid(indv):
    """Process offspring information when one of parents is haplodiploid"""
    exp_genosplito = set(indv[2].split(":")[0].split("/"))
    cont_genosplito = set(indv[3].split(":")[0].split("/"))
    if (len(exp_genosplito | cont_genosplito) == 2
            and exp_genosplito != cont_genosplito
            and (len(exp_genosplito) == 1
                 or len(cont_genosplito) == 1)):
        if len(exp_genosplito) == 1 and len(cont_genosplito) == 1:
            bothdiff = True
            exp_allele = list(exp_genosplito)[0]
            reverse = False
        else:
            bothdiff = False
            if (len(exp_genosplito) == 1
                    and ARGDICT["haplodiploid"] in ARGDICT["control_parent"]):
                exp_allele = list(exp_genosplito)[0]
                reverse = False
            elif (len(cont_genosplito) == 1
                  and ARGDICT["haplodiploid"] in ARGDICT["selected_parent"]):
                exp_allele = list(cont_genosplito)[0]
                reverse = True
            else:
                return(None)
        #### now let's get the offspring
        sample_scores = []
        for ns in range(2):
            sample_ix = list(set([call for call in indv[ns].split(":")[0].split("/")
                              if call == exp_allele]))
            if sample_ix:
                sample_ix = int(sample_ix[0])
                sample_score = float(indv[ns].split(":")[1].split(",")
                                 [sample_ix])/spt_vcfcov(indv[ns])
            else:
                sample_score = 0.0
            sample_scores.append(sample_score)
        high_scores = [i for i in sample_scores if i >= ARGDICT["mac"]]
        if (bothdiff or
                (not bothdiff
                 and len(high_scores) < len(sample_scores))):
            if reverse:
                sample_scores = [1-i for i in sample_scores]
            return(sample_scores)

def both_fixed(indv):
    """Process offspring information when both parents fixed"""
    exp_genosplito = set(indv[1].split(":")[0].split("/"))
    cont_genosplito = set(indv[2].split(":")[0].split("/"))
    if (len(exp_genosplito | cont_genosplito) == 2
            and exp_genosplito != cont_genosplito
            and (len(exp_genosplito) == 1
                 and len(cont_genosplito) == 1)):
        exp_allele = list(exp_genosplito)[0]
        sample_ix = list(set([call for call in indv[0].split(":")[0].split("/")
                              if call == exp_allele]))
        if sample_ix:
            sample_ix = int(sample_ix[0])
            sample_score = float(indv[0].split(":")[1].split(",")
                                 [sample_ix])/spt_vcfcov(indv[0])
        else:
            sample_score = 0.0
        return(sample_score)

def process_noparents(line, vcfline, cov):
    """Processes each line of a VCF file when parents are not specified"""
    # now let's parse it for each replicate:
    outdict = {}
    info = line[9:]
    for sample in zip(ARGDICT["selected_offspring"], ARGDICT["control_offspring"]):
        scores = []
        indv = [info[vcfline.index(sample[0])],
                info[vcfline.index(sample[1])]]
        if (len([call for call in indv if "." not in call]) == len(indv) and
                (cov[sample[0]]*ARGDICT["coverage_under"]
                 <= spt_vcfcov(indv[0]) <= cov[sample[0]]*ARGDICT["coverage_over"]
                    and cov[sample[1]]*ARGDICT["coverage_under"]
                    <= spt_vcfcov(indv[1]) <= cov[sample[1]]*ARGDICT["coverage_over"])):
            ########### NOW ONTO THE ACTUAL CALLS
            for call in indv:
                sample_score = float(call.split(":")[1].split(",")[0])/spt_vcfcov(call)
                scores.append(sample_score)
            high_scores = [i for i in scores if i >= ARGDICT["mac"]]
            low_scores = [i for i in scores if i <= 1-ARGDICT["mac"]]
            if (len(high_scores) < len(scores)
                    and len(low_scores) < len(scores)):
                outdict[sample[0]] = abs(scores[0] - scores[1])
                outdict[sample[1]] = 0
    return(outdict)

def process_infer(line, vcfline, cov):
    """Processes each line of a VCF file when only one parent is specified"""
    # now let's parse it for each replicate:
    outdict = {}
    info = line[9:]
    for sample in zip(ARGDICT["selected_offspring"],
                      ARGDICT["control_offspring"],
                      ARGDICT["major_parent"]):
        indv = [info[vcfline.index(sample[0])],
                info[vcfline.index(sample[1])],
                info[vcfline.index(sample[2])]]
        if (len([call for call in indv if "." not in call]) == len(indv) and
                (cov[sample[0]]*ARGDICT["coverage_under"]
                 <= spt_vcfcov(indv[0]) <= cov[sample[0]]*ARGDICT["coverage_over"]
                    and cov[sample[1]]*ARGDICT["coverage_under"]
                    <= spt_vcfcov(indv[1]) <= cov[sample[1]]*ARGDICT["coverage_over"]
                    and cov[sample[2]]*ARGDICT["coverage_under"]
                    <= spt_vcfcov(indv[2]) <= cov[sample[2]]*ARGDICT["coverage_over"])):
                ########### NOW ONTO THE ACTUAL CALLS
            sample_scores = inferred(indv)
            if sample_scores:
                outdict[sample[0]] = sample_scores[0]
                outdict[sample[1]] = sample_scores[1]
    return(outdict)

def process_hpd(line, vcfline, cov):
    """Processes each line of a VCF file when parents are not specified"""
    # now let's parse it for each replicate:
    outdict = {}
    info = line[9:]
    for sample in zip(ARGDICT["selected_offspring"],
                      ARGDICT["control_offspring"],
                      ARGDICT["selected_parent"],
                      ARGDICT["control_parent"]):
        indv = [info[vcfline.index(sample[0])],
                info[vcfline.index(sample[1])],
                info[vcfline.index(sample[2])],
                info[vcfline.index(sample[3])]]
        if (len([call for call in indv if "." not in call]) == len(indv) and
                (cov[sample[0]]*ARGDICT["coverage_under"]
                 <= spt_vcfcov(indv[0]) <= cov[sample[0]]*ARGDICT["coverage_over"]
                    and cov[sample[1]]*ARGDICT["coverage_under"]
                    <= spt_vcfcov(indv[1]) <= cov[sample[1]]*ARGDICT["coverage_over"]
                    and cov[sample[2]]*ARGDICT["coverage_under"]
                    <= spt_vcfcov(indv[2]) <= cov[sample[2]]*ARGDICT["coverage_over"]
                    and cov[sample[3]]*ARGDICT["coverage_under"]
                    <= spt_vcfcov(indv[3]) <= cov[sample[3]]*ARGDICT["coverage_over"])):
            ########### NOW ONTO THE ACTUAL CALLS
            sample_scores = haplodiploid(indv)
            if sample_scores:
                outdict[sample[0]] = sample_scores[0]
                outdict[sample[1]] = sample_scores[1]
    return(outdict)

def process_samples(line, vcfline, cov):
    """Processes each line of a VCF file"""
    # now let's parse it for each replicate:
    outdict = {}
    info = line[9:]
    for sample in zip(ARGDICT["selected_offspring"]+ARGDICT["control_offspring"],
                      ARGDICT["selected_parent"]*2,
                      ARGDICT["control_parent"]*2):
        indv = [info[vcfline.index(sample[0])],
                info[vcfline.index(sample[1])],
                info[vcfline.index(sample[2])]]
        if (len([call for call in indv if "." not in call]) == len(indv) and
                (cov[sample[0]]*ARGDICT["coverage_under"]
                 <= spt_vcfcov(indv[0]) <= cov[sample[0]]*ARGDICT["coverage_over"]
                 and cov[sample[1]]*ARGDICT["coverage_under"]
                 <= spt_vcfcov(indv[1]) <= cov[sample[1]]*ARGDICT["coverage_over"]
                 and cov[sample[2]]*ARGDICT["coverage_under"]
                 <= spt_vcfcov(indv[2]) <= cov[sample[2]]*ARGDICT["coverage_over"])):
                ########### NOW ONTO THE ACTUAL CALLS
            sample_score = both_fixed(indv)
            if sample_score is not None:
                outdict[sample[0]] = sample_score
    return(outdict)

def line_parser(line):
    """Parses VCF string to get mapping quality info"""
    linedict = {}
    line = line.split(";")
    for keyval in line:
        keyval = keyval.split("=")
        try:
            linedict[keyval[0]] = float(keyval[1])
        except ValueError:
            linedict[keyval[0]] = str(keyval[1])
    return(linedict)

def get_vcftuple():
    """Parses the VCF file, filters SNPs, and extracts relevant information"""
    print("PARSING VCF TO ANALYZE VARIANTS")
    chrom_tuple = scale_ends(ARGDICT["outdir1"] + "/chrom_file.txt")
    chrom_dict = chrom_tuple[0]
    chroms = chrom_tuple[1]
    number_total_snps = 0
    number_qc_snps = 0
    number_passed_snps = 0
    cov = {}
    with open(ARGDICT["outdir1"] + "/coverageinfo.txt", "r") as opencov:
        for line in opencov:
            line = (line.rstrip()).split("\t")
            strain = line[0]
            cov[strain] = float(line[1])
    masking = masker()
    openvcf = open(ARGDICT["vcf"], "r") # feeds the file line by line (see below)
    vcfdict = {}
    contigs = []
    vcfline = None
    outkeys = [i for i in ARGDICT["selected_offspring"]+ARGDICT["control_offspring"]]
    for line in openvcf: # let's go over each line
        line = line.replace("|", "/")
        line = line.rstrip().split("\t")
        chrom = line[0]
        vcfline = line[9:] if chrom == "#CHROM" else vcfline
        if (chrom in chroms
                and len(line[3]) == 1 and len(line[4]) == 1):
            pos = int(line[1])
            quals = line_parser(line[7])
            number_total_snps += 1
            if ("QD" in quals and "MQ" in quals and "SOR" in quals
                    and "MQRankSum" in quals and "ReadPosRankSum" in quals):
                if chrom not in contigs:
                    vcfdict[chrom] = {}
                    for binny in range(1, chrom_dict[chrom], ARGDICT["binsize"]):
                        vcfdict[chrom][binny] = {}
                    contigs.append(chrom)
                    current_bin = 1
                if pos >= current_bin + ARGDICT["binsize"]:
                    while pos >= current_bin + ARGDICT["binsize"]:
                        current_bin = current_bin + ARGDICT["binsize"]
                if (pos not in masking[chrom]
                        and quals["QD"] >= ARGDICT["qds"]
                        and quals["MQ"] >= ARGDICT["mps"]
                        and quals["SOR"] < ARGDICT["sor"]
                        and quals["MQRankSum"] >= ARGDICT["mqrs"]
                        and quals["ReadPosRankSum"] >= ARGDICT["rprs"]):
                    vcfdict[chrom][current_bin][pos] = {}
                    number_qc_snps += 1
                    if ("selected_parent" and "control_parent" in ARGDICT
                            and "haplodiploid" not in ARGDICT):
                        outdict = process_samples(line, vcfline, cov)
                    elif ("selected_parent" and "control_parent" in ARGDICT
                            and "haplodiploid" in ARGDICT):
                        outdict = process_hpd(line, vcfline, cov)
                    elif "major_parent" in ARGDICT:
                        outdict = process_infer(line, vcfline, cov)
                    else:
                        outdict = process_noparents(line, vcfline, cov)
                    if outdict:
                        number_passed_snps += 1
                        for sample in outdict:
                            vcfdict[chrom][current_bin][pos][sample] = outdict[sample]
    if "verbose" in ARGDICT:
        verbosy(vcfdict, contigs, outkeys)
    print("the total number of SNPs considered is %s"%(number_total_snps))
    print("the total number of SNPs that passed QC is %s"%(number_qc_snps))
    print("the total number of SNPs that passed QC and BSA cutoffs is %s"%(number_passed_snps))
    return(vcfdict, contigs, outkeys)

def process_segment(vcf_tuple, scaffy, beg):
    """Retrives allele counts within a genomic window"""
    relevant_keys = [i for i in vcf_tuple[0][scaffy] if i <= beg < i+ARGDICT["binsize"]
                     or beg <= i <= beg + ARGDICT["window"]
                     or i <= beg + ARGDICT["window"] < i+ARGDICT["binsize"]]
    segment_spls = {}
    for spl in vcf_tuple[2]:
        segment_spls[spl] = {}
        segment_spls[spl]["count"] = 0.0
        segment_spls[spl]["averages"] = 0.0
    for relevant_key in relevant_keys:
        section = sorted(vcf_tuple[0][scaffy][relevant_key])
        for snppy in section:
            if beg <= snppy <= beg + ARGDICT["window"]:
                snppy_spls = vcf_tuple[0][scaffy][relevant_key][snppy]
                for spl in snppy_spls:
                    spl_value = vcf_tuple[0][scaffy][relevant_key][snppy][spl]
                    segment_spls[spl]["count"] += 1
                    segment_spls[spl]["averages"] += spl_value
    return(segment_spls)

def slider(vcf_tuple):
    """Performs a sliding window analysis"""
    print("RUNNING SLIDING WINDOW ANALYSIS")
    shader = scale_dict(ARGDICT["outdir1"] + "/chrom_file.txt")[0]
    lennies = length_dict(ARGDICT["outdir1"] + "/chrom_file.txt")
    outdict = {}
    for spl in vcf_tuple[2]:
        outdict[spl] = {}
        outdict[spl]["val"] = []
        outdict[spl]["pos"] = []
        outdict[spl]["nvr"] = []
    for scaffy in vcf_tuple[1]:
        beg = 0
        medpos = mean([beg, beg+ARGDICT["window"]])+shader[scaffy]
        while (beg+ARGDICT["window"]) <= (lennies[scaffy] + ARGDICT["slide"]):
            segment_spls = process_segment(vcf_tuple, scaffy, beg)
            for spl in segment_spls:
                if segment_spls[spl]["count"] >= ARGDICT["min_allele"]:
                    outdict[spl]["val"].append(segment_spls[spl]["averages"]
                                               /segment_spls[spl]["count"])
                    outdict[spl]["pos"].append(medpos)
                    outdict[spl]["nvr"].append(int(segment_spls[spl]["count"]))
            beg = beg + ARGDICT["slide"]
            medpos = mean([beg, beg + ARGDICT["window"]]) + shader[scaffy]
    outdir2 = ARGDICT["outdir2"]
    if not os.path.isdir("%s"%ARGDICT["outdir2"]):
        subprocess.call("mkdir %s"%(outdir2), shell=True)
    if "verbose" in ARGDICT:
        for spl in vcf_tuple[2]:
            bsa_out = open(ARGDICT["outdir2"]+"/%s.txt"%(spl), "w")
            positions = outdict[spl]["pos"]
            for possa in positions:
                where = positions.index(possa)
                bsa_out.write(
                    "%s\t%s\t%s\n" %
                    (possa, outdict[spl]["val"][where],
                     outdict[spl]["nvr"][where]))
            bsa_out.close()
    return(outdict)

def shade_grid(axes, maxx):
    """Shades chromosomes in alternating gray and white"""
    shader_file = ARGDICT["outdir1"] + "/chrom_file.txt"
    shade_scaff_line = scale_dict(shader_file)
    shader = shade_scaff_line[0]
    scaffolds = shade_scaff_line[1]
    # here comes the shading of scaffolds/chroms
    old = None
    shade = False
    # here comes the shading
    for contig in scaffolds:
        val = shader[contig]
        if old and shade:
            axes.axvspan(old, val, color='0.87', alpha=0.5)
            shade = False
        else:
            if old != None:
                shade = True
        old = val
    # the last one
    if shade:
        axes.axvspan(old, maxx, color='0.87', alpha=0.5)
    axes.grid(True, axis="y", ls="dotted", color="0.85",lw=1.5)

def inducer(beg, end, number):
    """Spaces ticks on plot when spacing defined"""
    count = 0
    while number > 1e3:
        number = number/1e3
        count += 1
    ticklocator = arange(beg, end+number, number*math.pow(1e3, count))[1:]
    xminorlocator = None
    yminorlocator = None
    if ARGDICT["xminor"] > 0:
        xminorlocator = arange(beg, end+number, number*math.pow(1e3, count)/ARGDICT["xminor"])[1:]
    if ARGDICT["yminor"] > 0:
        yminorlocator = arange(beg, end+number, number*math.pow(1e3, count)/ARGDICT["yminor"])[1:]
    beglab = beg*math.pow(1e3, -count)
    endlab = end*math.pow(1e3, -count)
    ticklabels = arange(beglab, endlab+number, number)
    sigfigs = max([len(set(str(i))-set(["."])) for i in ticklabels])
    ticklabels = [round_sig(i, sigfigs) for i in ticklabels][1:]
    if ticklabels[-1] > 1e3:
        newbeg = beglab
        newend = endlab
        while ticklabels[-1] > 1e3:
            newbeg = newbeg/1e3
            newend = newend/1e3
            number = number/1e3
            count += 1
            ticklabels = arange(newbeg, newend+number, number)
    final_ticklabels = []
    for tickl in ticklabels:
        if int(tickl) == float(tickl):
            final_ticklabels.append(str(int(tickl)))
        else:
            final_ticklabels.append(str(tickl))
    return(ticklocator, final_ticklabels, xminorlocator, yminorlocator, count)

def process_coords(coord, number, count):
    """Processes coordinates for plotting"""
    if number < 1:
        number = "%f"%(number)
        number = number.rstrip("0")
        zeros = str(number).split(".")[1].count("0")
        number = float(number)
        newcoord = str(coord).split(".")
        newcoord[1] = newcoord[1][:zeros+1]
        newcoord = float(".".join(newcoord))
    else:
        newcoord = coord*math.pow(1e3, -count)
        newcoord = math.floor(newcoord)
        newcoord = newcoord*math.pow(1e3, count)
    return(newcoord, number)

def reducer(beg, end, nticks):
    """Helps determine tick spacing based on optimal number of ticks"""
    count = 0
    number = end - beg
    number = number/float(nticks)
    while number > 1e3:
        number = number/1e3
        count += 1
    number = float("%0.1g"%(number)) # math.floor(number) if number >= 1 else
    if number < 1:
        newbeg = process_coords(beg, number, count)[0]
        newend = process_coords(end, number, count)[0]
        number = process_coords(end, number, count)[1]
    else:
        newbeg = process_coords(beg, number, count)[0]
        newend = process_coords(end, number, count)[0]
    ticklocator = arange(newbeg, newend+number, number*math.pow(1e3, count))
    ticklocator = [i for i in ticklocator if beg < i < end]
    xminorlocator = None
    yminorlocator = None
    if ARGDICT["xminor"] > 0:
        xminorlocator = arange(newbeg, newend+number, number*math.pow(1e3, count)/ARGDICT["xminor"])
    if ARGDICT["yminor"] > 0:
        yminorlocator = arange(newbeg, newend+number, number*math.pow(1e3, count)/ARGDICT["yminor"])
    beglab = newbeg*math.pow(1e3, -count)
    endlab = newend*math.pow(1e3, -count)
    ticklabels = arange(beglab, endlab+number, number)
    if ticklabels[-1] > 1e3:
        newbeg = beglab
        newend = endlab
        while ticklabels[-1] > 1e3:
            newbeg = newbeg/1e3
            newend = newend/1e3
            number = number/1e3
            count += 1
            ticklabels = arange(newbeg, newend+number, number)
    ticklabels = [i for i in ticklabels if beg < int(i*math.pow(1e3, count)) < end]
    final_ticklabels = []
    for tickl in ticklabels:
        if int(tickl) == float(tickl):
            final_ticklabels.append(str(int(tickl)))
        else:
            final_ticklabels.append(str(tickl))
    return(ticklocator, final_ticklabels, xminorlocator, yminorlocator, count)

def ticks(axes, x_min, x_max, y_min, y_max):
    """Specifies ticks for the plot"""
    # x-axis
    if ARGDICT["xstep"]:
        scaled_ticks = arange(0, x_max+ARGDICT["xstep"], ARGDICT["xstep"])
        scaled_ticks = [i for i in scaled_ticks
                        if x_min <= i <= x_max+ARGDICT["xstep"]]
        tick_def = inducer(min(scaled_ticks), max(scaled_ticks), ARGDICT["xstep"])
    else:
        try:
            tick_def = reducer(x_min, x_max, ARGDICT["xticks"])
        except KeyError:
            error("YOU MUST PROVIDE EITHER TICK SPACING OR THE PREFERRED NUMBER OF TICKS")
    axes.set_xticks(tick_def[0])
    axes.set_xticklabels(tick_def[1], size=17.24)
    axes.set_xlabel("Genomic Position (%s)" %(PREFIX[tick_def[4]]), size=17.24)
    if ARGDICT["xminor"] > 0:
        axes.set_xticks(tick_def[2], minor=True)
    if ARGDICT["ystep"]:
        scaled_ticks = arange(-1.1, 1.1, ARGDICT["ystep"])
        scaled_ticks = [i for i in scaled_ticks if
                        y_min-ARGDICT["ystep"] <= i <= y_max+ARGDICT["ystep"]]
        tick_def = inducer(min(scaled_ticks),
                           max(scaled_ticks), ARGDICT["ystep"])
    else:
        try:
            tick_def = reducer(y_min, y_max, ARGDICT["yticks"])
        except KeyError:
            error("YOU MUST PROVIDE EITHER TICK SPACING OR THE PREFERRED NUMBER OF TICKS")
    axes.set_yticks(tick_def[0])
    if ARGDICT["yminor"] > 0:
        axes.set_yticks(tick_def[3], minor=True)
    axes.set_yticklabels(tick_def[1], size=17.24)
    axes.tick_params(axis="both", which="both", top=True, right=True)

def finish_plot(axes, indict, fig, commas):
    """Adds the option to zoom in on parts of plot"""
    all_spls = sorted(indict)
    if "nvr" in indict[all_spls[0]]:
        plot_file = ARGDICT["outdir3"] + "/BSA_sep_plot.pdf"
        axes.set_ylabel("Selected parent allele frequency", size=17.24)
    elif commas:
        plot_file = ARGDICT["outdir3"] + "/BSA_comb_plot.pdf"
        axes.set_ylabel("Difference in allele frequency", size=17.24)
    else:
        plot_file = ARGDICT["outdir3"] + "/BSA_average_plot.pdf"
        axes.set_ylabel("Difference in allele frequency", size=17.24)
    print("Saving plot to %s"%(plot_file))
    print("***********************")
    fig.savefig(plot_file, bbox_inches='tight')
    if "zoom_file" in ARGDICT:
        converta = scale_dict(ARGDICT["outdir1"] + "/chrom_file.txt")[0]
        zoomy = open(ARGDICT["zoom_file"])
        for line in zoomy:
            line = line.rstrip().split("\t")
            chrom = line[0]
            beg = int(line[1])+converta[chrom]
            end = int(line[2])+converta[chrom]
            y_min = 1
            y_max = 0
            for spl in all_spls:
                posses = [p for p in indict[spl]["pos"] if beg <= p <= end]
                where_beg = indict[spl]["pos"].index(posses[0])
                where_end = indict[spl]["pos"].index(posses[-1])
                values = indict[spl]["val"][where_beg:where_end]
                y_max = max(values) if max(values) > y_max else y_max
                y_min = min(values) if min(values) < y_min else y_min
            if "stat_cutoff" in ARGDICT:
                axes.axhline(y=ARGDICT["stat_cutoff"],
                             ls="dashed", lw=1, color=ARGDICT["sigcolor"])
                axes.axhline(y=-ARGDICT["stat_cutoff"],
                             ls="dashed", lw=1, color=ARGDICT["sigcolor"])
                y_max = ARGDICT["stat_cutoff"] if ARGDICT["stat_cutoff"] > y_max else y_max
                y_min = -ARGDICT["stat_cutoff"] if -ARGDICT["stat_cutoff"] < y_min else y_min
            new_extremes = (y_min-abs(y_max-y_min)*0.05, y_max+abs(y_max-y_min)*0.05)
            if "nvr" in indict[all_spls[0]]:
                y_min = 0
                y_max = 1
            else:
                y_min = new_extremes[0]
                y_max = new_extremes[1]
            ticks(axes, beg, end, y_min, y_max)
            if "nvr" not in indict[all_spls[0]]:
                axes.set_ylim(y_min, y_max)
            else:
                axes.set_ylim(-0.02, 1.02)
            axes.set_xlim(beg, end)
            new_plotfile = ARGDICT["outdir3"]+"/%s_%s_%s_%s.pdf"%(
                plot_file.split("/")[-1].split(".")[0], chrom, beg, end)
            print("Saving plot to %s"%(new_plotfile))
            print("***********************")
            fig.savefig(new_plotfile, bbox_inches='tight')

def plotter(indict):
    """Uses matplotlib to plot BSA scans"""
    ## FILES TO PLOT
    fig = plt.figure(figsize=(15, 5), dpi=1500)
    axes = plt.axes()
    print("***********************")
    print("PLOTTIN'")
    chrom_tuple = scale_ends(ARGDICT["outdir1"] + "/chrom_file.txt")
    if not os.path.isdir("%s"%ARGDICT["outdir3"]):
        subprocess.call("mkdir %s"%(ARGDICT["outdir3"]), shell=True)
    chrom_dict = chrom_tuple[0]
    chroms = chrom_tuple[1]
    x_min = 0
    x_max = chrom_dict[chroms[-1]]
    y_min = 1
    y_max = 0
    axes.axhline(y=0, ls="solid", lw=0.5, color="black")
    commas = None
    all_spls = sorted(indict)
    if ARGDICT["color"][0] != "custom":
        linespace = [0.01]+afill(0.01, 0.99, len(all_spls)-2)+[0.99]
        alterna = "plt.cm.%s(%s)"%(ARGDICT["color"], linespace)
        list_alt = list(eval(alterna))
        iter_alt = iter(list_alt)
    else:
        if "plot" not in ARGDICT:
            list_alt = ARGDICT["color"][1:3]
            iter_alt = iter(ARGDICT["color"][3:])
        else:
            list_alt = ARGDICT["color"][1:]
            iter_alt = iter(list_alt)
    for spl in all_spls:
        commas = True if "," in spl else False
        x_max = max(indict[spl]["pos"]) if max(indict[spl]["pos"]) > x_max else x_max
        x_min = min(indict[spl]["pos"]) if min(indict[spl]["pos"]) < x_min else x_min
        y_max = max(indict[spl]["val"]) if max(indict[spl]["val"]) > y_max else y_max
        y_min = min(indict[spl]["val"]) if min(indict[spl]["val"]) < y_min else y_min
        if "," in spl or spl == "average" or "perm" in indict[spl]:
            tsvet = next(iter_alt)
            linew = 1.5 if "," in spl or spl == "average" or "plot" in ARGDICT else 1
        else:
            tsvet = list_alt[1] if spl in ARGDICT["selected_offspring"] else list_alt[-2]
            linew = 1
        if "perm" in indict[spl]:
            y_max = indict[spl]["perm"] if indict[spl]["perm"] > y_max else y_max
            y_min = -indict[spl]["perm"] if -indict[spl]["perm"] < y_min else y_min
            axes.axhline(y=-indict[spl]["perm"], ls="dashed", lw=1, color=tsvet)
            axes.axhline(y=indict[spl]["perm"], ls="dashed", lw=1, color=tsvet)
        if ARGDICT["color"][0] != "custom":
            print("%s : %s" %(matplotlib.colors.to_hex(tsvet), spl))
        else:
            print("%s : %s" %(tsvet, spl))
        axes.plot(indict[spl]["pos"], indict[spl]["val"], lw=linew, color=tsvet)
    if "stat_cutoff" in ARGDICT:
        axes.axhline(y=ARGDICT["stat_cutoff"], ls="dashed", lw=1, color=ARGDICT["sigcolor"])
        axes.axhline(y=-ARGDICT["stat_cutoff"], ls="dashed", lw=1, color=ARGDICT["sigcolor"])
        y_max = ARGDICT["stat_cutoff"] if ARGDICT["stat_cutoff"] > y_max else y_max
        y_min = -ARGDICT["stat_cutoff"] if -ARGDICT["stat_cutoff"] < y_min else y_min
    # ACCESSORY STUFF TO MAKE IT PRETTY
    shade_grid(axes, x_max)
    new_extremes = (y_min-abs(y_max-y_min)*0.05, y_max+abs(y_max-y_min)*0.05)
    if "nvr" not in indict[all_spls[0]]:
        y_min = new_extremes[0]
        y_max = new_extremes[1]
    else:
        y_min = 0
        y_max = 1
    ticks(axes, x_min, x_max, y_min, y_max)
    # HOW FAR WILL THE PLOT GO?
    if "nvr" not in indict[all_spls[0]]:
        axes.set_ylim(y_min, y_max)
    else:
        axes.set_ylim(-0.02, 1.02)
    axes.set_xlim(x_min, x_max)
    finish_plot(axes, indict, fig, commas)

def plot_perm(filey):
    """Loads permutation information from file to plot"""
    cutoff = 0.0
    openfile = open(filey)
    for line in openfile:
        line = line.rstrip().split("\t")
        current = float(line[1])
        if current > cutoff:
            cutoff = current
    return(cutoff)

def plot_dict():
    """Loads sliding window information to plot"""
    outplot_dict = {}
    for filey in enumerate(ARGDICT["plot"]):
        outplot_dict["%s,%s"%(filey[0], filey[1])] = {}
        outplot_dict["%s,%s"%(filey[0], filey[1])]["pos"] = []
        outplot_dict["%s,%s"%(filey[0], filey[1])]["val"] = []
        if "permplot" in ARGDICT:
            outplot_dict["%s,%s"%(filey[0], filey[1])]["perm"] = plot_perm(
                ARGDICT["permplot"][filey[0]])
        openfile = open(filey[1])
        for line in openfile:
            line = line.rstrip().split("\t")
            outplot_dict["%s,%s"%(filey[0], filey[1])]["pos"].append(float(line[0]))
            outplot_dict["%s,%s"%(filey[0], filey[1])]["val"].append(float(line[1]))
    return(outplot_dict)

def lowest_highest(indict):
    """Defines highest and lowest limits for each chromosome"""
    beginnings = set()
    ends = set()
    for spl in indict:
        positions = indict[spl]["pos"]
        minpos = min(positions)
        maxpos = max(positions)
        beginnings.add(minpos)
        ends.add(maxpos)
    count = 0
    cur_beg = max(beginnings)
    cur_end = min(ends)
    while len(beginnings) != 1 or len(ends) != 1:
        beginnings = set([])
        ends = set([])
        for spl in indict:
            positions = [i for i in indict[spl]["pos"] if cur_beg <= i <= cur_end]
            minpos = min(positions)
            maxpos = max(positions)
            beginnings.add(minpos)
            ends.add(maxpos)
        cur_beg = max(beginnings)
        cur_end = min(ends)
        count += 1
        if count >= 100:
            error("NOT ENOUGH VALUES. CONSIDER RELAXING QC SETTINGS.")
    lowest = list(beginnings)[0]
    highest = list(ends)[0]
    return(lowest, highest)

def fill_loop(indict, outdict, spl, positions):
    """Loops over the internal part
     of sliding window positions and
     fills in missing positions"""
    for pos in positions:
        if pos != positions[0]:
            where = indict[spl]["pos"].index(pos)
            diff = indict[spl]["pos"][where] - indict[spl]["pos"][where-1]
            if diff == ARGDICT["slide"]:
                outdict[spl]["pos"].append(pos)
                outdict[spl]["val"].append(indict[spl]["val"][where])
                outdict[spl]["nvr"].append(indict[spl]["nvr"][where])
            if diff > ARGDICT["slide"]:
                nfills = int((diff/float(ARGDICT["slide"]))-1)
                outdict[spl]["pos"] = (outdict[spl]["pos"] +
                                       afill(indict[spl]["pos"][where-1],
                                             indict[spl]["pos"][where], nfills))
                outdict[spl]["val"] = (outdict[spl]["val"] +
                                       afill(indict[spl]["val"][where-1],
                                             indict[spl]["val"][where], nfills))
                outdict[spl]["nvr"] = outdict[spl]["nvr"] + [-1] * nfills
                outdict[spl]["pos"].append(pos)
                outdict[spl]["val"].append(indict[spl]["val"][where])
                outdict[spl]["nvr"].append(indict[spl]["nvr"][where])
        else:
            outdict[spl]["pos"].append(pos)
            outdict[spl]["val"].append(indict[spl]["val"][0])
            outdict[spl]["nvr"].append(indict[spl]["nvr"][0])
    return(outdict)

def fill_in(indict):
    """Makes sure every sample
     has the same number of values
     for permutations"""
    print("FILLING MISSING VALUES")
    outdict = {}
    chrom_tuple = scale_ends(ARGDICT["outdir1"] + "/chrom_file.txt")
    chrom_dict = chrom_tuple[0]
    chroms = chrom_tuple[1]
    for spl in indict:
        outdict[spl] = {}
        outdict[spl]["pos"] = []
        outdict[spl]["val"] = []
        outdict[spl]["nvr"] = []
    lo_hi = lowest_highest(indict)
    for spl in indict:
        positions = [pos for pos in indict[spl]["pos"]
                     if lo_hi[0] <= pos <= lo_hi[1]]
        # identify gaps in windows and fill them
        outdict = fill_loop(indict, outdict, spl, positions)
        # final fill
        diff = chrom_dict[chroms[-1]]-lo_hi[1]
        end_fills = 0
        if diff > ARGDICT["slide"]:
            fills = arange(lo_hi[1]+ARGDICT["slide"], chrom_dict[chroms[-1]], ARGDICT["slide"])
            end_fills = len(fills)
            outdict[spl]["pos"] = outdict[spl]["pos"]+fills
        diff = lo_hi[0]-(ARGDICT["slide"])
        beg_fills = 0
        if diff > ARGDICT["slide"]:
            fills = arange(0, ARGDICT["window"]/2,
                           ARGDICT["slide"])[1:]+arange(ARGDICT["window"]/2,
                                                        lo_hi[0], ARGDICT["slide"])
            beg_fills = len(fills)
            outdict[spl]["pos"] = fills + outdict[spl]["pos"]
        nfills = beg_fills+end_fills
        fill_val = []
        if nfills:
            fill_val = afill(indict[spl]["val"][-1], indict[spl]["val"][0], nfills)
        if fill_val:
            outdict[spl]["val"] = outdict[spl]["val"] + fill_val[:end_fills]
            outdict[spl]["val"] = fill_val[end_fills:] + outdict[spl]["val"]
            outdict[spl]["nvr"] = outdict[spl]["nvr"] + [-1] * end_fills
            outdict[spl]["nvr"] = [-1] * beg_fills + outdict[spl]["nvr"]
    for spl in outdict:
        bsa_out = open(ARGDICT["outdir2"]+"/%s_filled_in.txt"%(spl), "w")
        positions = outdict[spl]["pos"]
        for possa in positions:
            where = positions.index(possa)
            bsa_out.write(
                "%s\t%s\t%s\n" %
                (possa, outdict[spl]["val"][where], outdict[spl]["nvr"][where]))
        bsa_out.close()
    return(outdict)

def permute_shuffle(new_permute_dict):
    """Performs sliding permutations on replicates"""
    permuted_values = []
    # do i need this variable?
    ngrps = float(len(new_permute_dict))
    for grp in new_permute_dict:
        vals = new_permute_dict[grp]["val"]
        nvals = len(vals)
        where = random.randint(0, nvals)
        permvals = vals[where:] + vals[:where]
        if not permuted_values:
            permuted_values = permvals
        else:
            new_permuted_values = []
            for permval in zip(permvals, permuted_values):
                totes = sum(permval)
                new_permuted_values.append(totes)
            permuted_values = new_permuted_values
    permuted_values = [i/ngrps for i in permuted_values]
    tops = [abs(max(permuted_values)), abs(min(permuted_values))]
    topdist = max(tops)
    return(topdist)

def permute_process(new_permute_dict):
    """Performs series of sliding permutations"""
    perm_out = open(ARGDICT["outdir2"]+"/permutations.txt", "a")
    topdists = []
    for _perm in range(ARGDICT["perm"]):
        topdist = permute_shuffle(new_permute_dict)
        topdists.append(topdist)
    cutoff = (1-ARGDICT["sig"])*100
    critical_val = percentile(topdists, cutoff)
    outstring = ";".join(sorted(new_permute_dict))
    perm_out.write("%s\t%s\n"%(outstring, critical_val))
    perm_out.close()
    return(critical_val)

def unpermute(indict):
    """Average values among replicates"""
    noperm_dict = {}
    unpermuted_values = []
    for grp in indict:
        if not unpermuted_values:
            unpermuted_values = indict[grp]["val"]
            unpermuted_pos = indict[grp]["pos"]
        else:
            new_unpermuted_values = []
            for val in zip(indict[grp]["val"], unpermuted_values):
                totes = sum(val)
                new_unpermuted_values.append(totes)
            unpermuted_values = new_unpermuted_values
    unpermuted_values = [i/float(len(indict)) for i in unpermuted_values]
    noperm_dict["average"] = {}
    noperm_dict["average"]["val"] = unpermuted_values
    noperm_dict["average"]["pos"] = unpermuted_pos
    bsa_out = open(ARGDICT["outdir2"]+"/selected_average.txt", "w")
    positions = noperm_dict["average"]["pos"]
    for possa in positions:
        where = positions.index(possa)
        bsa_out.write("%s\t%s\n"%
                      (possa, noperm_dict["average"]["val"][where]))
    return(noperm_dict)

def combino(sel, unsel):
    """Combines specified selected and unselected samples"""
    outdict = {}
    outdict["pos"] = ARGDICT["master_dict"][sel]["pos"]
    outdict["val"] = []
    for val in zip(ARGDICT["master_dict"][sel]["val"], ARGDICT["master_dict"][unsel]["val"]):
        outdict["val"].append(val[0]-val[1])
    return(outdict)

# it permutes each replicate once and then gets the sum
def permute_setup():
    """Master function for permutation testing on BSA peaks"""
    perm_out = open(ARGDICT["outdir2"]+"/permutations.txt", "w")
    perm_out.close()
    dictlist = []
    comb_dict = {}
    for strains in zip(ARGDICT["selected_offspring"], ARGDICT["control_offspring"]):
        comb_dict["%s,%s"%(strains[0], strains[1])] = combino(strains[0], strains[1])
    dictlist.append(comb_dict)
    noperm_dict = unpermute(comb_dict)
    print("unpermuted min is %s"%(min(noperm_dict["average"]["val"])))
    print("unpermuted max is %s"%(max(noperm_dict["average"]["val"])))
    if ARGDICT["perm"] > 0:
        if "unpaired" in ARGDICT:
            print("COLLECTING DATA TO PERMUTE")
            new_comb_dict = {}
            unsels = list(permutations(ARGDICT["control_offspring"]))
            unsels = unsels[1:]
            random.shuffle(unsels)
            for unsel in unsels[:ARGDICT["combinations"]-1]:
                for strains in zip(ARGDICT["selected_offspring"], unsel):
                    new_comb_dict["%s,%s"%(strains[0],
                                           strains[1])] = combino(strains[0], strains[1])
                dictlist.append(new_comb_dict)
                new_comb_dict = {}
        print("RUNNING PERMUTATIONS")
        pool = multiprocessing.Pool(processes=ARGDICT["n_threads"])
        combocrit = pool.map(permute_process, dictlist)
        final_val = max(combocrit)
        # this is for plotting
        print("statistical cutoff is %s"%(final_val))
    else:
        final_val = None
    return(noperm_dict, comb_dict, final_val)

################################################################################################

if ("selected_offspring" in ARGDICT
        and "control_offspring" in ARGDICT):

    if ("control_parent" in ARGDICT
            and "selected_parent" in ARGDICT):

        ALL_GROUPS = [ARGDICT["selected_offspring"], ARGDICT["control_offspring"],
                      ARGDICT["control_parent"], ARGDICT["selected_parent"]]
        LENGTHS = set([len(ARGDICT["selected_offspring"]), len(ARGDICT["control_offspring"]),
                       len(ARGDICT["control_parent"]), len(ARGDICT["selected_parent"])])

    elif "major_parent" in ARGDICT:
        ALL_GROUPS = [ARGDICT["selected_offspring"], ARGDICT["control_offspring"],
                      ARGDICT["major_parent"]]
        LENGTHS = set([len(ARGDICT["selected_offspring"]), len(ARGDICT["control_offspring"]),
                       len(ARGDICT["major_parent"])])
    else:
        ALL_GROUPS = [ARGDICT["selected_offspring"], ARGDICT["control_offspring"]]
        LENGTHS = set([len(ARGDICT["selected_offspring"]), len(ARGDICT["control_offspring"])])
        if ARGDICT["perm"]:
            error("CANNOT RUN PERMUTATIONS IF PARENTAL DATA UNAVAILABLE")
            sys.exit()

    if len(LENGTHS) != 1:
        if len(LENGTHS) == 2 and 1 in LENGTHS:
            VALUE = int(list(LENGTHS - set([1]))[0])
            for group in ALL_GROUPS:
                if len(group) == 1:
                    for _ngroup in range(VALUE-1):
                        group.append(group[0])
        else:
            error("DIFFERENT NUMBERS OF INDIVIDUALS IN GROUPS. EXITING PROGRAM.")

    # finds average genome-wide read coverage for each strain/population in VCF file
    coverage()

    # goes over VCF and outputs allele frequencies to be used in sliding window analysis
    VCFTUPLE = get_vcftuple()

    # this outouts sliding windows
    FINAL_DICT = slider(VCFTUPLE)

    if not os.path.isdir("%s"%ARGDICT["outdir3"]):
        subprocess.call("mkdir %s"%(ARGDICT["outdir3"]), shell=True)

    NEW_FINAL_DICT = fill_in(FINAL_DICT)
    plotter(NEW_FINAL_DICT)
    ARGDICT["master_dict"] = NEW_FINAL_DICT
    PERM_RESULTS = permute_setup()
    AV_UNPERM_DICT = PERM_RESULTS[0]
    UNPERM_DICT = PERM_RESULTS[1]
    if PERM_RESULTS[2]:
        ARGDICT["stat_cutoff"] = PERM_RESULTS[2]
    plotter(AV_UNPERM_DICT)
    plotter(UNPERM_DICT)
elif "plot" in ARGDICT:
    PLOTTING_INFO = plot_dict()
    plotter(PLOTTING_INFO)

else:
    error("PLEASE EITHER PROVIDE PARENTAL/OFFSPRING OR PLOTTING INFORMATION")

#############################################################################################
