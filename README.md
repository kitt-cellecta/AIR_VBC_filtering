# Cellecta DriverMap AIR V2 - VBC Analysis <br />

The two scripts below both use the validation barcodes (VBCs) in the Cellecta DriverMap AIR to (1) filter clonotypes based on read counts and (2) estimate the number of template molecules (DNA or RNA) in the sample.

---

### Barcode Hopping and VBC filtering <br />

```bash
python3 run_filter_vbc.py DIRECTORY SAMPLE_NAME PERCENT_THRESHOLD
```
&nbsp;&nbsp;&nbsp;&nbsp; where: <br />
&nbsp;&nbsp;&nbsp;&nbsp; - DIRECTORY = directory of mixcr folder containing clones tsv files <br />
&nbsp;&nbsp;&nbsp;&nbsp; - SAMPLE_NAME = name of the current sample being analyzed <br />
&nbsp;&nbsp;&nbsp;&nbsp; - PERCENT_THRESHOLD = threshold value used during barcode hopping filtering; this value is the percentage; this value is multiplied to the most highly expressed read count value of a clonotype; VBCs with less value than this multiplication product are deemed as barcode hoppers and are filtered out from the rest of the analysis <br />
e.g.
```bash
python3 run_filter_vbc.py ~/Desktop/test_T/ D_T-CDR3-VBC-2ug 5
```

---

### Template Molecule Estimation

```bash
python3 quantify_templates.py DIRECTORY SAMPLE_NAME
```
&nbsp;&nbsp;&nbsp;&nbsp; where: <br />
&nbsp;&nbsp;&nbsp;&nbsp; - DIRECTORY = directory of mixcr folder containing clones tsv files <br />
&nbsp;&nbsp;&nbsp;&nbsp; - SAMPLE_NAME = name of the current sample being analyzed <br />
e.g.
```bash
python3 quantify_templates.py ~/Desktop/test_T/ D_T-CDR3-VBC-2ug
```
