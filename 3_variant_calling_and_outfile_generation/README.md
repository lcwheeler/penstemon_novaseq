## Genotyping and variant calling


### In testing phase
This directory is still in testing and optimization phase. Check back soon ;)


samtools faidx $reference genome



## check the formatting for the other READMEs


For filtering. Thoughts:
(1) It's OK to call reference when Qual is low.
(2) When Qual is low for a variant site, either change to missing or call as reference. I think it's fine to call as reference in these cases. So probably you'd want to flag the variant as lowqual, and then in the filtering step, change all lowqual variants to reference.

*both of these changes might require changes to the filtering script*

(3) Compare the output of the vcf to the gvcf. Want to make sure that there isn't something lost in the conversion step.