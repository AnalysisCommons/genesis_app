#!/bin/bash
# CHARGE_SeqMeta2 0.0.2
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

main() {

    echo "Value of phenofile: '$phenofile'"
    echo "Value of snpinfofile: '$snpinfofile'"
    echo "Value of genotypefile: '$genotypefile'"
    echo "Value of outputfilename: '$outputfilename'"

    echo "Value of kinshipmatrix: '$kinshipmatrix'"
    echo "Value of pheno_id: '$pheno_id'"
    echo "Value of snpNames: '$snpNames"
    echo "Value of nsmatch: '$nsmatch'"

    echo "Value of buffer: '$buffer'"
    echo "Value of genefile: '$genefile'"
    echo "Value of snp_filter: '$snp_filter'"
    echo "Value of gene_filter: '$gene_filter'"
    echo "Value of top_maf: '$top_maf'"
    covariate_list=$(echo ${covariate_list} | sed 's/ //g') 
    echo "Value of covariate_list: '$covariate_list'"
    

    #set -x
    #  dx download "$phenofile" -o phenofile
    #  dx download "$snpinfofile" -o snpinfofile
    
    # COMMENTS: MB 12/4/15
    # Loop calls process() to run genotype files in parallel
    # 'genotypefile' used to be an array of Rdata files, but is now just the one gds file (see dxapp.json)
    # No changes made yet to process() to reflect this because still working locally to get the DNANexusAnalysisPipeline_test.R file running
    # Of course, this affect postprocess() too, but I haven't changed that yet either.
    # Maybe instead of looping through genotype files we loop through chr number



    dx download "$phenofile" -o phenofile &
    dx download "$genotypefile" -o genotypefile &
    dx download "$snpinfofile" -o snpinfofile &
    dx download "$kinshipmatrix" -o kinshipmatrix &
    dx download "$genefile" -o genefile &
    # install R
    # add R to path 
    ls
    
    echo "INSTALLING GENESIS"
    make & 
    export PATH=/opt/R/bin/:${PATH}
    export MKL_NUM_THREADS=1
    wait
    ls
##   echo "FIXING PHENO"
## Fixing pheno file 
##    head -n1 phenofile > fixed_pheno
##    grep -v NA phenofile | \
##         tail -n +2 | sort -k1,1 | \
##         awk -F, 'BEGIN{OFS=",";prev=""} $1!=prev {prev = $1;print $0}' >> fixed_pheno
##    sed -i -e 's/"f"/"F"/' -e 's/"m"/"M"/' fixed_pheno
##    mv fixed_pheno phenofile
    sudo chmod o+rw /tmp
    # wait if debug 
    if [ ${debug} -ne 0 ]
    then
       echo "DEBUG is on sleeping for ${debug}h"
       echo "Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" snpinfofile genotypefile results $kinshipmatrix $pheno_id  $nsmatch $buffer genefile $snp_filter $gene_filter $top_maf  $test_requested $burden_test $min_mac"
       sleep ${debug}h
    fi
    wait
    echo "Checking phenofile" 
    if [ -e phenofile ] 
    then
       head -n1 phenofile
    else
       echo "The phenofile is not ready"
    fi
 
    echo "Running code"
    Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" snpinfofile genotypefile results kinshipmatrix $pheno_id  $nsmatch $buffer genefile $snp_filter $gene_filter $top_maf  $test_requested $burden_test $min_mac
    echo "Finished running code"
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    dx mv ${results} ${outputfilename}.csv.gz
}
