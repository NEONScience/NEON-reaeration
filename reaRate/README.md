NEON Reaeration Rate
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- ****** Description ****** -->

Contains functions to format data and calculate gas-exchange rates of
tracer gas (SF6) and oxygen.

<!-- ****** Usage ****** -->

## Usage

The functions in this package have the following purposes: format
downloaded data, calculate SF6 loss rates, calculate travel time, and
normalize loss rate to Schmidt 600 number. See help files for individual
functions for details and the exampleScript.R file in the testingAndDev
folder of this repo for a working example. The general flow of using
this package is:

1.  Load data into R environment (the examples use
    neonUtilities::loadByProduct)
2.  Run def.format.reaeration, which returns a data frame with data from
    many tables merged by station and date (usually 4 rows per NaCl or
    NaBr experiment).
3.  Run gas.loss.rate.plot, which outputs many QAQC plots into specified
    directory and returns a dataframe with SF6 loss rates using raw gas
    data, outlier removed gas data, outlier removed gas data corrected
    for plateau salt concentration, and outlier removed gas data
    corrected with background subtracted plateu salt concentration.
4.  Optional, run bkgd.salt.conc.plot to see plots of background salt
    concentrations over time and flow and relationships between sensor
    data and grab sample data.
5.  Run def.calc.reaeration, which returns a list of two dataframes, the
    input dataframe of data for up to 4 stations per site per date and
    an output dataframe appended with loss rate, travel time, SF6
    reaeration rate, O2 gas transfer velocity, and Schmidt number 600
    for a given site and date.

<!-- ****** Acknowledgements ****** -->

## Credits & Acknowledgements

<!-- HTML tags to produce image, resize, add hyperlink. -->
<!-- ONLY WORKS WITH HTML or GITHUB documents -->

<a href="http://www.neonscience.org/">
<img src="logo.png" width="300px" /> </a>

<!-- Acknowledgements text -->

The National Ecological Observatory Network is a project solely funded
by the National Science Foundation and managed under cooperative
agreement by Battelle. Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author(s)
and do not necessarily reflect the views of the National Science
Foundation.

<!-- ****** License ****** -->

## License

GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

<!-- ****** Disclaimer ****** -->

## Disclaimer

*Information and documents contained within this package are available
as-is. Codes or documents, or their use, may not be supported or
maintained under any program or service and may not be compatible with
data currently available from the NEON Data Portal.*
