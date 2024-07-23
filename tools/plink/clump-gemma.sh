filt_file=${1?Please provide PLINK file root as first parameter}

report=${2?Please provide GWAS report file as second parameter}

rep=${report##*/}

outbase="cl_${rep%.assoc.*}"

rangefile=${4:-"dummy.range"}

# Clumping parameters

# Top peak pvalue thresh
#p1=${4:-"1.2e-4"}
p1=${3:-"1e-6"}

# Minimum peak pvalue thresh
#p2="8e-4"
p2=${5:-"0.001"}

# Minimum r2 for inclusion in the peak
r2=${6:-"0.5"}

# Maximum gap (flanking region around the ranges in the rangefile), kb
gap_kb=${7:-"20"}

# Maximum distance for inclusion in the peak
kb=${8:-"250"}


# Concise clumping report, with ranges
plink --allow-no-sex --bfile $filt_file   --clump $report  \
     --clump-field p_wald --clump-snp-field rs  \
	 --keep-allele-order  \
	 --clump-p1 $p1   \
	 --clump-p2 $p2  \
	 --clump-kb $kb --clump-r2 $r2  \
	 --out $outbase  \
	 --clump-range $rangefile   --clump-range-border $gap_kb

rm ${outbase}.nosex

# For verbose report (using ranges)
plink --allow-no-sex --bfile $filt_file  --clump $report  \
    --clump-field p_wald --clump-snp-field rs    --keep-allele-order \
	--clump-p1 $p1  \
	 --clump-p2 $p2 \
	 --clump-kb $kb --clump-r2 $r2  \
	 --out ${outbase}.verbose  \
	 --clump-verbose  \
	 --clump-range $rangefile  --clump-range-border $gap_kb 

rm ${outbase}.verbose.nosex
