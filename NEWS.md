# GencoDymo2 1.0.4

## Bug fixes
* Added missing global variable declaration to fix dplyr update (PR #04, thanks @DavisVaughan).


# GencoDymo2 1.0.3

## Bug Fixes
- Fixed issue [#2](https://github.com/monahton/GencoDymo2/issues/2): Corrected handling of gene_type argument in the `compare_release()` function. Now the argument is gene_biotype

# GencoDymo2 1.0.2

## Bug Fixes
- Removed hard dependency on `BSgenome.Hsapiens.UCSC.hg38`; now listed in `Suggests`.
- Updated multiple functions to handle genome input more robustly using `requireNamespace()` pattern.
- Ensured CRAN compatibility for functions using `getSeq()`.

## Improvements
- Cleaned up examples to avoid running genome-heavy code during CRAN checks.

# GencoDymo2 1.0.1

* Initial stable release of GencoDymo2
* Provides functions to extract, compare, and analyze GENCODE annotations.
* Supports generation of splice site motif FASTA files.
* A modified remake of GencoDymo package

# GencoDymo2 1.0.0.9000

* Initial development release of GencoDymo2
* A modified remake of GencoDymo package
