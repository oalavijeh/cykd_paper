# Create the full functional variant list
# Make a two column TSV, the first column being the original variant ID, and the second column the new variant ID (all same chromosome, in ascending order)
# Use the full functional variant list to extract only the variants of interest from the PGEN files
# Merge these much smaller PGEN files into one file
# Create a new fileset with --make-pgen --update <file> using the TSV from step 2 as the input to update.
# Use this new fileset with the updated IDs as the input to the AVT workflow, along with the updated variant list
# To phrase differently, you'll need to:
#
# extract the variants of interest from the plink files
# merge these files together
# rename all the variants so they are on one chromosome
# rerun the workflow with the new plink data + the new variant list (both with renamed variants)

module load bio/PLINK/2.00-devel-20200409-x86_64
for i in {1..22} X; do plink2 --make-bed --pfile chr${i}_masked_post_protein --out chr${i}_masked_post_protein; done

module load bio/PLINK/1.9b_4.1-x86_64
ls chr*.pgen | sed 's/.pgen//g' > to_merge.txt
plink --bfile chr1_masked_post_protein --merge-list to_merge.txt --make-bed --out merged_post_protein
plink --bfile merged_post_protein --update-name variant_input_sorted_for_relabel.txt --make-bed --out merged_post_protein_relabelled --allow-no-sex

module load bio/PLINK/2.00-devel-20200409-x86_64
