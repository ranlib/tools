#
# v1.0.5:
#
# variant_interpretation
# add rrs gene condition
# comment column -> warning column, fill with comment depending on %region > threshold
#
# vcf->dataframe
# filter information list as ; separated string
# new input: bed file
#
# v1.0.6
# add check that there is annotation for a variant at all, ie there is an ANN key
# add option to filter variants if no PASS in filter vcf volumn
# add verbose flag, message on stdout can be captured
#
# v1.0.7
# change vcf_to_pandas_dataframe: choice of annotation from snpEff
# choose 1st annotation if not modifier
# if modifier variant, loop over all annotations and choose closest
#
# v1.0.8
# add variant QC step for variants at end of variant_interpretation:
# add QC failed to warning column
# do not remove all large deletions, take all deletions for now
#
# v1.0.9
# add new column mdl_LIMSfinal to lab report
# copy mdl severity into mdl_LIMSfinal,
# update that severity bases on variant QC and regions coverage QC
#
# v1.1.0
# revamp filter_genes section: Gene_Name can be &-separated list of genes
# for large deletions
# move this section to beginning, so genes with no antimicrobial get
# antimicrobial set by using chemicals associated with region
# remove Warning column,
# set 'Insufficient Coverage' for mdl_LIMSfinal severity if 'Breadth_of_coverage_QC' == FAIL
# except for variants with mdl severity R
# go back to old 2.2.1 but with feature_ablation added
# still check for upstream_gene_variant for is_notsynonymous
# change to csv report
#
# v1.2.0
# rename function 3_2 to 3_2_2 and update function:
# 3.2.2: for non-synonymous variants introduce distinction: large deletions vs rest
#
variant_interpretation:
	docker build --no-cache -t dbest/variant_interpretation:v1.2.0 -f Dockerfile.variant_interpretation .
	#docker push dbest/variant_interpretation:v1.0.9 # push new version
	#docker rmi dbest/variant_interpretation:v1.0.7 # remove old version


#
# v1.0.0
# initial version
#
# v1.0.1
# remove verbose flag, introduce logging
# set severity to WT for variants that contain failed QC in Warning column
# remove bed file, not needed
# take chemial and gene/chemical information from hard coded dictionaries
#
# v1.0.2
# filter out genes in input lab report that are not in the header of
# the lims report, remove duplicate terms in cells of lims report,
# turn lims report from tsv into csv file
#
lims_report:
	docker build --no-cache -t dbest/lims_report:v1.0.2 -f Dockerfile.lims_report .
	#docker push dbest/lims_report:v1.0.1 # push new version
	#docker rmi dbest/lims_report:v1.0.0 # remove old version

#
# v1.0.0
# initial version
# docker version of get_lineage
# == CDC lineage parser but with command line options and use of vcf file
#
lineage:
	docker build -t dbest/lineage:v1.0.0 -f Dockerfile.lineage .

#
# v1.0.0
# initial version of fast_lineage_caller
#
fast_lineage_caller:
	docker build -t dbest/fast_lineage_caller:v1.0.0 -f Dockerfile.fast_lineage_caller .
