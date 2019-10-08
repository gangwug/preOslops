#! /usr/local/bin/julia
############################################################
#introduction page for julia language https://docs.julialang.org/en/stable/
#parse the input parameters
using ArgParse     #using the module for command-line argument parsing, see examples at 'https://github.com/carlobaldassi/ArgParse.jl/tree/master/examples'
s = ArgParseSettings()
@add_arg_table s begin
    "--cycdir", "-c"
        help = "the complete path of CYCLOPS code. eg. /TOOLS/code/JuliaCYCLOPS2a"
	required = true
    "--infile", "-i"
        help = "the input expression dataset, csv file (including path)"
        required = true
    "--seedfile", "-s"
        help = "the seed gene list, csv file (including path)"
        required = true
    "--outdir", "-o"
        help = "the directory store ouput csv files"
        required = true
    "--Frac_Var"
	help = "Set Number of Dimensions of SVD to maintain this fraction of variance"
	arg_type = Float64
	default = 0.99
    "--DFrac_Var"
	help = "Set Number of Dimensions of SVD to so that incremetal fraction of variance of var is at least this much"
	arg_type = Float64
	default = 0.01
    "--Seed_MinCV"
	help = "Set the minimal CV for filtering the genes in the seed list, genes with CV below this value will be removed"
	arg_type = Float64
	default = 0.14
    "--Seed_MaxCV"
	help = "Set the maximal CV for filtering the genes in the seed list, genes with CV above this value will be removed"
	arg_type = Float64
	default = 0.85
    "--Seed_Blunt"
	help = "Set the blunt number for removing outliers"
	arg_type = Float64
	default = .975
    "--Seed_MinMean"
	help = "Set the mimimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed"
	arg_type = Float64
	default = -100000000000000.00
    "--Out_Symbol"
	help = "The symbol used in the output file"
end
#the actual argument parsing step, parsed_args is a dictionary type
parsed_args = parse_args(ARGS, s)
############################################################
addprocs(5)
using StatsBase
using MultivariateStats
using Distributions
#change to the directory contains CYCLOPS code and import required modules
cd(string(parsed_args["cycdir"]))
using CYCLOPS_2a_AutoEncoderModule
using CYCLOPS_2a_PreNPostprocessModule
using CYCLOPS_2a_MCA
using CYCLOPS_2a_MultiCoreModule_Smooth
using CYCLOPS_2a_Seed
############################################################
Seed_MinCV				= parsed_args["Seed_MinCV"]
Seed_MaxCV				= parsed_args["Seed_MaxCV"]
Seed_Blunt				= parsed_args["Seed_Blunt"]
Seed_MinMean				= parsed_args["Seed_MinMean"]
Frac_Var                                = parsed_args["Frac_Var"]               # Set Number of Dimensions of SVD to maintain this fraction of variance
DFrac_Var                               = parsed_args["DFrac_Var"]              # Set Number of Dimensions of SVD to so that incremetal fraction of variance of var is at least this much
###########################################################
#read the seed gene list
seedfile=string(parsed_args["seedfile"])
fullnonseed_data_BHTC=readcsv(seedfile)
bhtc_seeds=fullnonseed_data_BHTC[2:end ,1]

#get the eigen gene expression profile
fullnonseed_data_merge=readcsv(string(parsed_args["infile"]))
fullnonseed_data2=hcat(fullnonseed_data_merge[:,1], fullnonseed_data_merge[:,1], fullnonseed_data_merge)       ## the default get seed function assumes first column=probe, 2nd symbol, 3rd entrez (or just text)
alldata_samples=fullnonseed_data2[1,4:end]
seed_symbols_bhtc, seed_data_bhtc 					= getseed(fullnonseed_data2,bhtc_seeds,Seed_MaxCV,Seed_MinCV,Seed_MinMean,Seed_Blunt)
seed_data_bhtc 							        = dispersion!(seed_data_bhtc)
outs_bhtc, norm_seed_data_bhtc, eigen_sig_fraction, eigen_sum_fraction	= GetEigenGenes(seed_data_bhtc,Frac_Var,DFrac_Var,300)

#output the dispersed seed gene expression and the eigen expression matrix
seed_out = hcat(seed_symbols_bhtc, seed_data_bhtc)
seed_out = vcat(hcat("geneName", alldata_samples), seed_out)
writecsv(string(parsed_args["outdir"], "/", parsed_args["Out_Symbol"], "_SeedgeneExp.csv"), seed_out)
eigen_out = hcat( map(string, fill("eigen_", outs_bhtc), 1:outs_bhtc, fill("_", outs_bhtc), round(eigen_sig_fraction[1:outs_bhtc], 4) ), norm_seed_data_bhtc )
eigen_out = vcat(hcat("eigenName", alldata_samples), eigen_out)
writecsv(string(parsed_args["outdir"], "/", parsed_args["Out_Symbol"], "_EigengeneExp.csv"), eigen_out)
###########################################################
#write the file parameters out
outdir=string(parsed_args["outdir"])
parafile=string(parsed_args["Out_Symbol"], "_para.txt")
f = open(string(outdir, "/", parafile), "a")
for (key, val) in parsed_args 
    write(f, "$key : $val \n")
end
write(f, "\nThe important intermediate parameters in ordering with each eigen cluster is listed below:\n")
close(f)
