Usage:
python3 run_filter_vbc.py <dir> <sample_name> <percent_threshold> <br />
e.g <br />
python3 run_filter_vbc.py ~/Desktop/test_T/ D_T-CDR3-VBC-2ug 5 <br />
<br />
where: <br />
<dir> = directory of mixcr folder containing clones tsv files <br />
<sample_name> = name of the current sample being analyzed <br />
<percent_threshold> = threshold value used during barcode hopping filtering; <br />
	this value is the percentage; this value is multiplied to the most highly <br />
expressed read count value of a clonotype; VBCs with less value than this <br />
multiplication product are deemed as barcode hoppers and are filtered out <br />
from the rest of the analysis <br />
