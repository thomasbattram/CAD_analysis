#### put alleles_new on github inputs :) 
#### remove the data that is called here within the r scripts
#### add a bit at the end of each section to save the image and load again in next script


main_dir=$1
metab=$2
gen=$3
cohort_stats=$4

# main_dir="PATH_TO_DIR"
# metab="METABOLITE_DATA_FILE"
# gen="GENETIC_DATA_FILE"
# cohort_stats="COHORT_DATA_FILE"

cd "${main_dir}"

Rscript R/Masterfile.R "${metab}" "${gen}" "${cohort_stats}"
