#! /usr/local/bin/julia
############################################################
#introduction page for julia language https://docs.julialang.org/en/stable/
#parse the input parameters
using ArgParse     #using the module for command-line argument parsing, see examples at 'https://github.com/carlobaldassi/ArgParse.jl/tree/master/examples'
s = ArgParseSettings()
@add_arg_table s begin
    "--cycdir", "-c"
        help = "the complete path of CYCLOPS code (eg. /TOOLS/code/JuliaCYCLOPS2a)"
	required = true
    "--infile", "-i"
        help = "the input expression dataset, csv file (including path)"
        required = true
    "--normseed", "-n"
        help = "the complete file name (including path) of seed expression data outputed from 'runCYCLOPS_Eigen.jl' "
        required = true
    "--eigenfile", "-e"
        help = "the complete file name (including path) of eigen expression data outputed from 'runCYCLOPS_Eigen.jl' "
        required = true    
    "--oscopecluster", "-p"
        help = "the complete file name (including path) of eigen expression data outputed from 'runCYCLOPS_Eigen.jl' "
        required = true    
    "--outdir", "-o"
        help = "the directory store CYCLOPS' ouput csv files"
        required = true
    "--Seed_Random"
	help = "the random seed number used for CYCLOPS"
	arg_type = Int
	default  = 12345
    "--N_best"
	help = "Number of random initial conditions to try for each optimization"
	arg_type = Int
	default  = 40
    "--N_background"
	help = "Number of background runs for global background distribution"
	arg_type = Int
	default  = 40
    "--Seed_MinMean"
	help = "Set the mimimal mean expression cutoff for seed list, genes with mean expression level below this value will be removed"
	arg_type = Float64
	default = -100000000000000.00
    "--Out_Symbol"
	help = "The symbol used in the output file to differentiate multiple runs for the same input file"
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
N_best                                  = parsed_args["N_best"]                 # Number of random initial conditions to try for each optimization
total_background_num                    = parsed_args["N_background"]           # Number of background runs for global background ref
Seed_MinMean				= parsed_args["Seed_MinMean"]
Seed_Random                             = parsed_args["Seed_Random"]
###########################################################
#write the file parameters out
outdir = string(parsed_args["outdir"])
parafile = string(parsed_args["Out_Symbol"], "_para.txt")
f = open(string(outdir, "/", parafile), "a")
write(f, string("CYCLOPS StartTime:", strftime(time()), "\nThe input file parameters for running CYCLOPS is as below: \n"))
for (key, val) in parsed_args 
    write(f, "$key : $val \n")
end
write(f, "\nThe important intermediate parameters in analyzing each data file is listed below:\n")
close(f)
###########################################################
#read in the expression data
fullnonseed_data_merge=readcsv(string(parsed_args["infile"]))
fullnonseed_data2=hcat(fullnonseed_data_merge[:,1], fullnonseed_data_merge[:,1], fullnonseed_data_merge) ## the default get seed function assumes first column=probe, 2nd symbol, 3rd entrez (or just text)
alldata_samples=fullnonseed_data2[1,4:end]

#read in the data output from 'runCYCLOPS_Eigen.jl'
seed_data_bhtc = readcsv(string(parsed_args["normseed"]))
seed_symbols_bhtc = seed_data_bhtc[2:end, 1]
seed_data_bhtc = float64(seed_data_bhtc[2:end, 2:end])
norm_seed_data_bhtc = readcsv(string(parsed_args["eigenfile"]))
norm_seed_data_bhtc = float64(norm_seed_data_bhtc[2:end,2:end])
###########################################################
srand(Seed_Random)
oscope_eigen_group = readcsv(string(parsed_args["oscopecluster"]))
for i in 2:( size(oscope_eigen_group, 1) )
    tindex = split(oscope_eigen_group[i, 3], "|")  
    tindex = map(x -> parse(Int, x), tindex)  
    norm_seed_data_oscope = norm_seed_data_bhtc[tindex,:]
    outs_bhtc_oscope = length(tindex)
    estimated_phaselist_bhtc,bestnet_bhtc,global_var_metrics_bhtc 	= CYCLOPS_Order(outs_bhtc_oscope,norm_seed_data_oscope,N_best)
    estimated_phaselist_bhtc						= mod(estimated_phaselist_bhtc + 2*pi,2*pi)
    #do not change the 'norm_seed_data_bhtc' in the 'smoothness_measures' function, since it uses the first eigengene expression value to get a linear order
    global_smooth_metrics_bhtc						= smoothness_measures(seed_data_bhtc,norm_seed_data_bhtc,estimated_phaselist_bhtc)
    global_metrics_bhtc	     = global_var_metrics_bhtc
    pvals = 1
    ##note: the pvalue here is not completed same as the originally design for CYCLOPS, because CYCLOPS select the top eigengenes, so it is an equal comparison here
    ##after Oscope selection, the eigengenes used is not the top eigengenes, so the comparison here is not complete equal; this need to contact with Ron for detail discussion
    pvals = multicore_backgroundstatistics_global_eigen(seed_data_bhtc,outs_bhtc_oscope,N_best,total_background_num,global_metrics_bhtc)
    ############################################################
    #get the primary cosinor analysis results
    cosinor_skin_bhtc=Compile_MultiCore_Cosinor_Statistics(fullnonseed_data2, estimated_phaselist_bhtc, 4, 48)
    ##adjust the phase with E-box genes
    eboxgenes=["DBP", "HLF", "TEF", "PER1", "NR1D2", "CIART", "C1orf51", "PER3"]
    eboxphases_bhtc=cosinor_skin_bhtc[[findin(cosinor_skin_bhtc[:,2],eboxgenes)],:]
    bhtc_criteria=((eboxphases_bhtc[:,4].<.1) & ( eboxphases_bhtc[:,11].>1.25) & ( eboxphases_bhtc[:,9].>Seed_MinMean))
    eboxphases_bhtc=eboxphases_bhtc[findin(bhtc_criteria,true),:]
    ebox_genes = eboxphases_bhtc[:,2] 
    eboxphases_bhtc=float(eboxphases_bhtc[:,6])
    ############################################################
    function Circular_Mean(phases::Array{Float64,1})
      sinterm=sum(sin(phases))
      costerm=sum(cos(phases))
      atan2(sinterm,costerm)
    end
    ############################################################
    if (convert(Bool, length(eboxphases_bhtc)) ) 
	  phaseTag=["EboxAdjusted"]
	  eboxTag=join(ebox_genes, ",")
	  eboxphase_bhtc=Circular_Mean(eboxphases_bhtc)
	  estimated_phaselist_bhtc_final=mod(estimated_phaselist_bhtc .- eboxphase_bhtc+(pi),2*pi)  ##ebox genes in mouse lung peak on average CT 11.5
    else
	  phaseTag=["RawPhase"]
	  eboxTag="NaN"
	  estimated_phaselist_bhtc_final=estimated_phaselist_bhtc
    end
    # suppose that in the human study, there are more samples collected in the day time than the night time, which may be not always true in the real application
    BHTC_Nday=length(findin((estimated_phaselist_bhtc_final .>pi),true))
    BHTC_Nnight=length(findin((estimated_phaselist_bhtc_final .<pi),true))
    if (BHTC_Nday<BHTC_Nnight)
	estimated_phaselist_bhtc_final=mod(2*pi-estimated_phaselist_bhtc_final,2*pi)
    end
    ############################################################
    # Transcript Tracings
    cosinor_skin_bhtc=Compile_MultiCore_Cosinor_Statistics(fullnonseed_data2,estimated_phaselist_bhtc_final,4,48)
    out_pha=alldata_samples
    cnumber=length(out_pha) 
    out_pha=vcat(out_pha, reshape(estimated_phaselist_bhtc,1,cnumber), reshape(estimated_phaselist_bhtc_final,1,cnumber))
    out_pha=hcat(reshape(Any["SampleID" "Original_Phase" "EboxAdjusted_Phase"], 3, 1), out_pha)  
    #write out the cosinor analysis results
    writecsv(string(parsed_args["outdir"], "/", parsed_args["Out_Symbol"], oscope_eigen_group[i, 1],  ".ToRunCYCLOPS_PhaseCYCLOPS.csv"), out_pha)
    writecsv(string(parsed_args["outdir"], "/", parsed_args["Out_Symbol"], oscope_eigen_group[i, 1],  ".ToRunCYCLOPS_CosinorCYCLOPS.csv"), cosinor_skin_bhtc)
    #println(string("\nSmooth:", out_paras["global_smooth"], "\nPvalue:", out_paras["pvals"], "\n"))
    ############################################################
    #transfer the  parameters to 'out_paras' and store the important intermediate values into it
    out_paras                     = Dict()
    out_paras["cluster_name"]     = oscope_eigen_group[i,1]
    out_paras["eigen_name"]       = oscope_eigen_group[i,2]
    out_paras["global_smooth"]    = global_smooth_metrics_bhtc[1]
    out_paras["pvals"]            = pvals
    # write the parameters stored in out_paras out
    f =  open(string(outdir, "/", parafile), "a")
    for (key, val) in out_paras 
    	write(f, "$key : $val  \n")
    end
    write(f, string("\nThe output phase is: ", string(phaseTag), "\nThe genes used to adjust phase are: ", string(eboxTag), "\n") )
    close(f)
end
############################################################
# write the time of fnishing the analysis
f =  open(string(outdir, "/", parafile), "a")
write(f, string("CYCLOPS StopTime:", strftime(time()),"\n\n\n"))
close(f)
