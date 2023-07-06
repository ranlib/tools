# bgzip myVariants.vcf; bcftools index myVariants.vcf.gz
# Plotting the mutation distribution for 'NC_000962.3'
from pysam import VariantFile
import seaborn as sns
filename = "ERR552797_DR_loci_raw_annotation.vcf.gz"
vcf_in = VariantFile(filename)

x = []
for rec in vcf_in.fetch("NC_000962.3"):
    x.append(rec.pos)

sns.distplot(x)
