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
# v1.2.1
# a) changes to severity in 1.1, 2.1, 3.1
# b) include has-deletion condition in region_coverage_qc
# c) rename mdl_2_2 and looker_2_2 to mdl_2_1 and looker_2_1 to make naming more consistent
# d) remove feature_ablation condition in region_coverage_QC, not needed anymore since
#    has_deletion condition used now
# e) change term has_large_deletion to has_deletion since we take all
#    DELLY deletion which can be also smaller than 50bp
#    this need to be improved, use also mutect2 deletions
# f) better get_deletions_in_region script: now get deletion interval
#    and check for overlaps with array of regions of interest
#
#
# v1.3.0
# a) split large deletion across several genes into several rows, one row per gene
# b) revamp get deletions in region
# c) add # bp deletion takes in region + gene length of lab report as debug option
#
# v1.3.1
# a) create looker and mdl_prelim columns beforehand with defaults,
#    just like the coverage columns,
#    code just fills the new columns
#    those columns are the output produced
#
# v1.4.0
# 2.2.2: change CDS to AA,
# 2.2.1: use distance,
# 3.2.2: split into 3.2.2.1/2/3, special treatment of gyrA and gyrB, use AA.pos
# section 1: remove genes Rv0678, mmpL5, mmpS5 from gene_list_1 and create gene_list_4, special treatment of gene_list_4 in section 1
#
variant_interpretation:
	docker build --no-cache -t dbest/variant_interpretation:v1.4.0 -f Dockerfile.variant_interpretation .
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
# v1.0.3
# remove 'No mutations detected' if there are any mutations
#
# v1.0.4
# reset index after filtering for reportable genes
# check for empty input lab report, output empty lims report in this case
#
lims_report:
	docker build --no-cache -t dbest/lims_report:v1.0.4 -f Dockerfile.lims_report .
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
