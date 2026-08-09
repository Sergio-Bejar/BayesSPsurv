## Resubmission

This package was archived on 2023-06-14 because issues were not corrected in time.
All issues have now been addressed:

- Removed `CXX_STD = CXX11` from Makevars and Makevars.win (C++11 is now the default).
- Fixed `Rcpp` listed in Imports but not imported: added `importFrom(Rcpp, sourceCpp)`.
- Removed erroneous `export(.onUnload)` from NAMESPACE.
- Removed `src-i386/` directory containing compiled binary objects.
- Replaced invalid URLs (HTTP 403) with `\doi{}` references.
- Replaced `length(class(object))` with `inherits()` for safe class checking.
- Updated RoxygenNote to 7.3.2.

There is a new maintainer: Sergio Bejar <sergio.bejar@cide.edu>
(previous maintainer: Nicolas Schmidt <nschmidt@cienciassociales.edu.uy>).

## Test environments
* local macOS, R 4.x
* R-hub (Linux, Windows)
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes
