## Resubmission

This package was archived on 2023-06-14 because issues were not corrected in time.
All issues have now been addressed:

- Removed `CXX_STD = CXX11` from Makevars and Makevars.win (C++11 is now the default).
- Fixed `Rcpp` listed in Imports but not imported: added `importFrom(Rcpp, sourceCpp)`.
- Removed erroneous `export(.onUnload)` from NAMESPACE.
- Removed `src-i386/` directory containing compiled binary objects.
- Replaced invalid URLs (HTTP 403) with `\doi{}` references.
- Replaced `length(class(object))` with `inherits()` for safe class checking.
- Fixed code/documentation mismatches in `form` defaults and `parameter` defaults for summary methods.
- Fixed `bsps_map()` to convert column names to numeric before `countrycode` gwn lookup.
- Added documentation for new `dem_survival` dataset.
- Replaced deprecated `@docType package` with `"_PACKAGE"` sentinel.
- Regenerated all .Rd files with roxygen2 8.1.0.

There is a new maintainer: Sergio Bejar <sergio.bejar@cide.edu>
(previous maintainer: Nicolas Schmidt <nschmidt@cienciassociales.edu.uy>).

## Test environments
* local macOS (x86_64-apple-darwin17.0), R 4.2.1

## R CMD check results

0 errors | 0 warnings | 2 notes

* NOTE: New submission. Package was archived on CRAN on 2023-06-14.
* NOTE: One URL in README.md (https://doi.org/10.1177/0022002714528006) returns
  HTTP 403. This is a valid DOI; SAGE Journals blocks automated URL checkers.
