### Usage <br />
<br />
<br />
## Barcode Hopping and VBC filtering <br />
<br />
python3 run_filter_vbc.py \<directory\> \<sample\_name\> \<percent_threshold\> <br />
e.g <br />
python3 run_filter_vbc.py ~/Desktop/test_T/ D_T-CDR3-VBC-2ug 5 <br />
<br />
where: <br />
\<directory\> = directory of mixcr folder containing clones tsv files <br />
\<sample\_name\> = name of the current sample being analyzed <br />
\<percent\_threshold\> = threshold value used during barcode hopping filtering; <br />
	this value is the percentage; this value is multiplied to the most highly <br />
expressed read count value of a clonotype; VBCs with less value than this <br />
multiplication product are deemed as barcode hoppers and are filtered out <br />
from the rest of the analysis <br />
<br />
<br />
## Template Molecule Estimation
<br />
python3 quantify_templates.py \<directory\> \<sample\_name\> <br />
e.g <br />
python3 quantify_templates.py ~/Desktop/test_T/ D_T-CDR3-VBC-2ug <br />
<br />
where: <br />
\<directory\> = directory of mixcr folder containing clones tsv files <br />
\<sample\_name\> = name of the current sample being analyzed <br />