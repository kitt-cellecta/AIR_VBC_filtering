# Usage:
# python3 run_filter_vbc.py <DIRECTORY> <sample_name> <percent_threshold>
# e.g 
# python3 run_filter_vbc.py ~/Desktop/test_T/ D_T-CDR3-VBC-2ug 5
# 
# where:
# <DIRECTORY> = directory of mixcr folder containing clones tsv files
# <sample_name> = name of the current sample being analyzed
# <percent_threshold> = threshold value used during barcode hopping filtering;
# 	this value is the percentage; this value is multiplied to the most highly 
#	expressed read count value of a clonotype; VBCs with less value than this
#	multiplication product are deemed as barcode hoppers and are filtered out
#	from the rest of the analysis
