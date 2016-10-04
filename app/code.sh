#!/bin/bash
# CHARGE_SeqMeta2 0.0.2

main() {

    echo "Value of phenofile: '$phenofile'"
    echo "Value of snpinfofile: '$snpinfofile'"
    echo "Value of genotypefile: '$genotypefile'"
    echo "Value of outputfilename: '$outputfilename'"

    echo "Value of kinshipmatrix: '$kinshipmatrix'"
    echo "Value of pheno_id: '$pheno_id'"
    echo "Value of snpNames: '$snpNames"

    echo "Value of buffer: '$buffer'"
    echo "Value of genefile: '$genefile'"
    echo "Value of snp_filter: '$snp_filter'"
    echo "Value of gene_filter: '$gene_filter'"
    echo "Value of top_maf: '$top_maf'"

    echo "Value of min_mac: '$min_mac'"

    echo "Value of test_type: '$test_type'"
    echo "Value of test_stat: '$test_stat'"
    echo "Value of weights: '$weights'"


    covariate_list=$(echo ${covariate_list} | sed 's/ //g') 
    echo "Value of covariate_list: '$covariate_list'"
    


    dx download "$phenofile" -o phenofile &
    dx download "$genotypefile" -o genotypefile &
    dx download "$snpinfofile" -o snpinfofile &
    dx download "$kinshipmatrix"   &

    # genefile
    if [[ "$genefile" != "" ]] ; then
	echo 'downloading genefile'
	
	dx download "$genefile" -o genefile &
	ingenefile="genefile"
	#genefile = genefile
    else
	ingenefile="NO_GENE_REGION_FILE"
    fi

    # install R
    # add R to path 
    
    echo "INSTALLING GENESIS"
    make >> /dev/null & 
    export PATH=/opt/R/bin/:${PATH}
    export MKL_NUM_THREADS=1
    wait

    echo "\nKINSHIP"
    kinshipmatrix_filename=$( dx describe --name "$kinshipmatrix" )
    
 

    sudo chmod o+rw /tmp
    # wait if debug 
    if [ ${debug} -ne 0 ]
    then
       echo "DEBUG is on sleeping for ${debug}h"
       echo "Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" snpinfofile genotypefile results $kinshipmatrix_filename $pheno_id  $buffer $ingenefile $snp_filter $gene_filter $top_maf $test_requested $burden_test $min_mac $weights"
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
 
    echo "Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" snpinfofile genotypefile results $kinshipmatrix_filename $pheno_id  $buffer $ingenefile $snp_filter $gene_filter $top_maf $test_requested $burden_test $min_mac $weights"
    echo "Running code"
    Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" snpinfofile genotypefile results $kinshipmatrix_filename $pheno_id  $buffer $ingenefile $snp_filter $gene_filter $top_maf  $test_stat $test_type $min_mac $weights
    echo "Finished running code"
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    dx mv ${results} ${outputfilename}.csv.gz
}
