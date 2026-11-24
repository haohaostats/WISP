# GitHub release preparation

- MIT license: confirmed by the authors.
- Repository root: this WISP folder, not the manuscript or research project root.
- Authors: Hao Chen and Clara Grazian; maintainer email hche2036@uni.sydney.edu.au.
- Initial version: 0.1.0.
- Public API: wisp(), wmst_pseudo(), wisp_ipcw().
- Help: installed Rd pages and README; no PDF manual.
- Cases: installed scripts load breastCancer and tandmob from optional icensBKL; no real records are redistributed.
- Code style: no comments or roxygen blocks in any R code, including examples and tests.
- Continuous integration: R CMD check --no-manual on Windows, Linux and macOS, to run after upload.
- Local build products and check evidence: dist/ (ignored by Git and package build).

GitHub destination: haohaostats/WISP. Repository URLs are included in DESCRIPTION, README and citation metadata. Do not invent a paper DOI or publication status. Run CI after upload and tag v0.1.0 only when those checks pass. CRAN submission is a separate future task.
