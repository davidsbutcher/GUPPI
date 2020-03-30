Get UniProt Protein Info (GUPPI)
================

Process TDPortal top-down reports (tdReports) by retrieving information
from the UniProt webservice and filtering by a selectable FDR value.

## Input

The processing of tdReports is carried out by the `guppi()` function.
Arguments are as follows:

### Mandatory arguments

  - `filedir` Add the path to the folder containing tdReports to be
    processed. This directory can be at any level above the tdReports in
    the directory structure.

  - `filename` Full name of the tdReport file or files, including
    extension. Must match the name of file exactly.

  - `taxon_number` UniProt taxon number to use for retrieving UniProt
    data. Should match the taxon number in the tdReport file.

  - `outputdir` Directory to place output files. Must be an existing
    folder.

### Optional arguments

  - `fdr` False detection rate to use for filtering results. Defaults to
    1%.

  - `make_dashboard` Boolean value (TRUE or FALSE) which controls
    whether an HTML report is generated.

  - `use_PB` Boolean value (TRUE or FALSE). If set to true,
    `RPushbullet` will be used to send a Pushbullet notification to your
    device when the analysis is finished. See `?RPushbullet::pbSetup`
    for more info.

## Analysis of tdReports

A connection is established to the SQLite database in the TD Report
using `RSQLite`. The “main” output includes all protein isoform
accession numbers with the lowest Q value from among all hits for each
isoform and the name of the data file from which the lowest Q value hit
was obtained. All isoforms with Q values which are missing or greater
than the cutoff value are deleted. Output is also generated which
contains all hits for all isoforms that are above the FDR cutoff
(`/allproteinhits`) and lowest Q-value hits sorted by data file
(`/proteinsbydatafile`).

UniProt data is added as detailed in the `Flat Files` section, with the
exception that average and monoisotopic masses are taken directly from
the tdReport file.

Minimum Q value from among all hits, average and monoisotopic masses,
and data file for lowest Q value hit are obtained for all proteoforms.
Proteoforms whose Q values are above the FDR cutoff are deleted.
Proteoforms whose corresponding protein isoform entry is above the FDR
cutoff are also deleted. UniProt info for every unique proteoform record
number is copied from corresponding protein entries to avoid wasting
time by querying the UniProt database again.

## Output

Output files are saved to the output directory (`outputdir`). Files are
timestamped with the time the script was initialized or share the same
name as the input file.

  - Protein results are saved to
    `outputdir/YYYYMMDD_hhmmss_protein_results.xlsx`.
  - Proteoform results are saved to
    `outputdir/YYYYMMDD_hhmmss_proteoform_results.xlsx`.
  - Output xlsx files contain sheets corresponding to each input
    tdReport file and a final summary sheet containing all input file
    names and counts of cytosolic, membrane, periplasmic, and
    unannotated proteins.
  - Lists of all protein hits including UniProt accession number, Q
    value, data file, result type (i.e. tight absolute mass, find
    unexpected modifications, or biomarker) are saved to
    `outputdir/protein_results_allhits` and share names with input
    files.
  - Workspace images (.Rdata file containing all R objects) from the end
    of the analysis are saved to `outputdir/workspace_images`.
  - If `make_dashboard` is set to TRUE, an HTML report is saved to
    `outputdir/report`.
