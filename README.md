# Cellecta DriverMap AIR V2 - VBC Analysis <br />

The two scripts below both use the validator barcodes (VBCs) in the Cellecta DriverMap AIR to (1) filter clonotypes based on read counts and (2) estimate the number of template molecules (DNA or RNA) in the sample. These two scripts are used sequentially and both operate directly on the MiXCR output folder. Template molecule estimation requires the results of the barcode hopping/VBC filtering steps.

---

### Barcode Hopping and VBC filtering <br />

```bash
python3 run_filter_vbc.py $DIRECTORY $SAMPLE_NAME --mode $MODE --peak-selection $PEAK_SELECTION
```
&nbsp;&nbsp;&nbsp;&nbsp; where: <br />
&nbsp;&nbsp;&nbsp;&nbsp; - DIRECTORY = directory of mixcr folder containing clones tsv files <br />
&nbsp;&nbsp;&nbsp;&nbsp; - SAMPLE_NAME = name of the current sample being analyzed (match with the MiXCR output prefix) <br />
&nbsp;&nbsp;&nbsp;&nbsp; - MODE = 'bulk' or 'single_cell' (Default: 'bulk') <br />
&nbsp;&nbsp;&nbsp;&nbsp; - PEAK_SELECTION = `by_distance` or `by_prominence` for KDE distributions with more than two peaks (Default: `by_distance`, which reproduces the original pipeline behavior) <br />
e.g.
```bash
python3 run_filter_vbc.py ~/Desktop/test_T/ test_T --mode bulk --peak-selection by_distance
```

`by_prominence` evaluates only adjacent peak pairs and selects the pair with
the strongest combination of local KDE prominence and clonotype constituent
support. It avoids selecting a sparse far-right tail simply because it creates
a favorable position-based split. The existing valley-depth check is applied
after either strategy selects a peak pair.

---

### Template Molecule Estimation

```bash
python3 run_quantify.py $DIRECTORY $SAMPLE_NAME --mode $MODE
```
&nbsp;&nbsp;&nbsp;&nbsp; where: <br />
&nbsp;&nbsp;&nbsp;&nbsp; - DIRECTORY = directory of mixcr folder containing clones tsv files <br />
&nbsp;&nbsp;&nbsp;&nbsp; - SAMPLE_NAME = name of the current sample being analyzed (match with the MiXCR output prefix) <br />
&nbsp;&nbsp;&nbsp;&nbsp; - MODE = 'bulk' or 'single_cell' (Default: 'bulk') <br />
e.g.
```bash
python3 run_quantify.py ~/Desktop/test_T/ test_T --mode bulk
```
