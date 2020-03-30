Get UniProt Protein Info (GUPPI)
================

Process TDPortal top-down reports (tdReports) by retrieving information
from the UniProt webservice and filtering by a selectable FDR value.

## Installation

Install from Github with:

``` r
remotes::install_github("davidsbutcher/GUPPI")
```

## Input

The processing of tdReports is carried out by the `guppi()` function. An
example of running the function:

``` r
guppi(
   "C:/Users/David Butcher/TDReports",
   c(
      "20200420_Excellent_TDReport_01.tdReport",
      "20200420_Excellent_TDReport_02.tdReport"
   ),
   83333,
   "C:/Users/David Butcher/Documents/guppi_output",
   fdr = 0.01,
   make_dashboard = TRUE,
   use_PB = FALSE
)
```

Arguments to the `guppi` function are as follows:

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
using `RSQLite`. All protein-level (or isoform) and proteform-level IDs
and other relevant data for each ID are extracted. The taxon number is
checked against files in the package directory to see if a corresponding
UniProt taxon database has already been downloaded. If not, the UniProt
web service is queried for all UniProt accession numbers in the taxon
using the package `UniProt.ws`. Protein name, organism, organism taxon
ID, protein sequence, protein function, subcellular location, and any
associated GO IDs are returned. Note that some of these values may not
be found and come back as empty or NA.

The UniProt taxon database is used to add information for all IDs
extracted from the tdReport. GO terms are obtained for all GO IDs using
the `GO.db` package and terms corresponding to subcellular locations are
saved in column “GO\_subcellular\_locations”. Average and monoisotopic
masses are determined from the protein sequence (intact sequence for
isoform-level IDs and proteoform sequence for proteoform-level IDs)
using the `Peptides` package.

Minimum Q value from among all hits, average and monoisotopic masses,
and data file for lowest Q value hit are obtained for all proteoforms.
Proteoforms whose Q values are above the FDR cutoff are deleted.
Proteoforms whose corresponding protein isoform entry is above the FDR
cutoff are also deleted, as in the TDViewer software.

## Output

The “main” output includes Q values, observed precursor masses, data
files, subcellular locations from the GO database and a variety of other
parameters for all isoform and proteoform IDs. All isoform IDs and
proteoform IDs with Q values which are missing or greater than the
cutoff value (`fdr`) are deleted.

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
