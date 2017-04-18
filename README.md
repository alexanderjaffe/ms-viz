# ms-viz
MS2 visualization and analysis tool to be used in conjunction with [mzStats](https://github.com/alexcritschristoph/mzStats/).

## Tutorial

### Download and configure

Download the ms-viz folder using the "Clone or Download button" on the top right.

In your terminal, navigate to the ms-viz-master folder and create a "data" directory.

```mkdir data```

Now, you're ready to run the pipeline.

### Process spectra

```
usage: process_spectra.py [-h] -i INPUT -c COMPOUND_TABLE -s SPECTRA_FILE

Converts a set of mzxml to json format for d3 viz and mgf files for gnps
querying.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to directory of mzxml.
  -c COMPOUND_TABLE, --compound_table COMPOUND_TABLE
                        Path to filtered compound table.
  -s SPECTRA_FILE, --spectra_file SPECTRA_FILE
                        Path to spectra mapping file.
```

For example: ```test```

### Query GNPS

```
usage: query_GNPS.py [-h] -i INPUT

Gets chemical info from GNPS for a set of spectra.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to directory of mgf.
```

### Format results

```
usage: format_results.py [-h] -i INPUT [-min MINIMUM_MZ] [-max MAXIMUM_MZ]

Reformats GNPS and compound table files for visualization.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to filtered compound table.
  -min MINIMUM_MZ, --minimum_mz MINIMUM_MZ
                        Minimum mz to keep.
  -max MAXIMUM_MZ, --maximum_mz MAXIMUM_MZ
                        Maximum mz to keep.
```

### Start the web visualization

In your terminal, type the command:

``` python -m SimpleHTTPServer 8001```

This starts up a local web-like server on which we will use the visualization.

Next, navigate to localhost:8001 on any web browser.
