#!/bin/bash

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
    echo "Value of varaggfile: '$varaggfile'"
    echo "Value of snp_filter: '$snp_filter'"
    echo "Value of gene_filter: '$gene_filter'"
    echo "Value of top_maf: '$top_maf'"

    echo "Value of min_mac: '$min_mac'"

    echo "Value of test_type: '$test_type'"
    echo "Value of test_stat: '$test_stat'"
    echo "Value of weights: '$weights'"


    covariate_list=$(echo ${covariate_list} | sed 's/ //g') 
    echo "Value of covariate_list: '$covariate_list'"
    echo "value of conditional SNP: '$conditional'"    
    echo "Value of user_cores: '$user_cores'" 
    echo "Value of het_vars: '$het_vars'"


    dx download "$phenofile" -o phenofile &
    dx download "$genotypefile" -o genotypefile &
    dx download "$kinshipmatrix"   &

    # snpinfofile
    if [[ "$snpinfofile" != "" ]] ; then
	echo 'downloading snpinfofile'
	
	dx download "$snpinfofile" -o snpinfofile &
	insnpinfofile="snpinfofile"
    else
	insnpinfofile="NO_SNPINFO_FILE"
    fi

    # genefile
    if [[ "$genefile" != "" ]] ; then
	echo 'downloading genefile'
	
	dx download "$genefile" -o genefile &
	ingenefile="genefile"
    else
	ingenefile="NO_GENE_REGION_FILE"
    fi


    # varaggfile
    if [[ "$varaggfile" != "" ]] ; then
	echo 'downloading varaggfile'
	
	dx download "$varaggfile" -o varaggfile &
	invaraggfile="varaggfile"
    else
	invaraggfile="NO_VAR_AGG_FILE"
    fi

    # install R
    # add R to path 
    
    echo "INSTALLING GENESIS"
    make >> /dev/null & 
    export PATH=/opt/R/bin/:${PATH}
    export MKL_NUM_THREADS=1
    wait

    echo "KINSHIP"
    kinshipmatrix_filename=$( dx describe --name "$kinshipmatrix" )
    
 

    sudo chmod o+rw /tmp
    # wait if debug 
    if [ ${debug} -ne 0 ]
    then
       echo "DEBUG is on sleeping for ${debug}h"
       echo "Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" $insnpinfofile genotypefile results $kinshipmatrix_filename $pheno_id  $buffer $ingenefile $invaraggfile $snp_filter $gene_filter $top_maf $test_requested $burden_test $min_mac $weights $user_cores $het_vars"
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
 
    echo "Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" $insnpinfofile genotypefile results $kinshipmatrix_filename $pheno_id  $buffer $ingenefile $invaraggfile $snp_filter $gene_filter $top_maf $test_stat $test_type $min_mac $weights $conditional $het_vars"
    echo "Running code"
    Rscript genesis.R phenofile $outcome_name $outcome_type \"$covariate_list\" $insnpinfofile genotypefile results $kinshipmatrix_filename $pheno_id  $buffer $ingenefile $invaraggfile $snp_filter $gene_filter $top_maf  $test_stat $test_type $min_mac $weights $conditional $user_cores $het_vars
    echo "Finished running code"
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    dx mv ${results} ${outputfilename}.csv.gz
}
