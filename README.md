# Cellecta DriverMap AIR V2 - VBC Analysis <br />

The two scripts below both use the validator barcodes (VBCs) in the Cellecta DriverMap AIR to (1) filter clonotypes based on read counts and (2) estimate the number of template molecules (DNA or RNA) in the sample. These two scripts are used sequentially and both operate directly on the MiXCR output folder. Template molecule estimation requires the results of the barcode hopping/VBC filtering steps.

---

### Complete VBC pipeline

Use the combined runner to perform filtering and then template-molecule
quantification for one sample:

```bash
bash run_vbc_pipeline.sh $DIRECTORY $SAMPLE_NAME --mode $MODE --peak-selection $PEAK_SELECTION
```

For example:

```bash
bash run_vbc_pipeline.sh ~/Desktop/test_T/ test_T --mode bulk --peak-selection by_prominence
```

The runner accepts `bulk`, `single_cell`, and `bead` modes and the
`by_distance`, `by_prominence`, and `manual` peak-selection modes. It stops if
filtering fails. For a manual rerun that creates a new filtered table, it
automatically refreshes template quantification instead of retaining an older
quantification report. Set `PYTHON_BIN` if the pipeline requires a particular
Python interpreter, for example `PYTHON_BIN=/path/to/python bash ...`.

---

### Barcode Hopping and VBC filtering <br />

```bash
python3 run_filter_vbc.py $DIRECTORY $SAMPLE_NAME --mode $MODE --peak-selection $PEAK_SELECTION
```
&nbsp;&nbsp;&nbsp;&nbsp; where: <br />
&nbsp;&nbsp;&nbsp;&nbsp; - DIRECTORY = directory of mixcr folder containing clones tsv files <br />
&nbsp;&nbsp;&nbsp;&nbsp; - SAMPLE_NAME = name of the current sample being analyzed (match with the MiXCR output prefix) <br />
&nbsp;&nbsp;&nbsp;&nbsp; - MODE = 'bulk', 'single_cell', or 'bead' (Default: 'bulk') <br />
&nbsp;&nbsp;&nbsp;&nbsp; - PEAK_SELECTION = `by_distance`, `by_prominence`, or `manual` (Default: `by_distance`, which reproduces the original pipeline behavior) <br />
e.g.
```bash
python3 run_filter_vbc.py ~/Desktop/test_T/ test_T --mode bulk --peak-selection by_distance
```

`by_prominence` evaluates only adjacent peak pairs and selects the pair with
the strongest combination of local KDE prominence and clonotype constituent
support. It avoids selecting a sparse far-right tail simply because it creates
a favorable position-based split. The existing valley-depth check is applied
after either strategy selects a peak pair.

### Manual valley selection

After a completed automatic run, the sample directory contains
`peak_manual_select.csv`, with eight VBC rows. Set active rows to `keep` or a
labelled KDE valley such as `valley_2`, then rerun with
`--peak-selection manual`. The manual run reads the original completed-run
manifest, archives automatic VBC output under `original_vbc_analysis/`, and
uses selected valleys as final filters without applying a second round of
post-selection threshold adjustments. In 4-VBC modes, rows 5–8 are ignored.

---

### Template Molecule Estimation

```bash
python3 run_quantify.py $DIRECTORY $SAMPLE_NAME --mode $MODE [--force]
```
&nbsp;&nbsp;&nbsp;&nbsp; where: <br />
&nbsp;&nbsp;&nbsp;&nbsp; - DIRECTORY = directory of mixcr folder containing clones tsv files <br />
&nbsp;&nbsp;&nbsp;&nbsp; - SAMPLE_NAME = name of the current sample being analyzed (match with the MiXCR output prefix) <br />
&nbsp;&nbsp;&nbsp;&nbsp; - MODE = 'bulk', 'single_cell', or 'bead' (Default: 'bulk') <br />
&nbsp;&nbsp;&nbsp;&nbsp; - `--force` = regenerate template estimates when a prior template-estimation report exists <br />
e.g.
```bash
python3 run_quantify.py ~/Desktop/test_T/ test_T --mode bulk
```
