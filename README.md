

# Meta-analysis of changes in phenotypic means and variation in urban populations

This repository contains materials used for a project investigating changes in phenotypic means and variation in urban bird populations.

---

Pablo Capilla-Lasheras, Megan J. Thompson, Alfredo Sánchez-Tójar, Yacob Haddou, Claire J. Branston, Denis Réale, Anne Charmantier, Davide M. Dominoni. **A global meta-analysis reveals higher variation in breeding phenology in urban birds than in their non-urban neighbours**. *bioRxiv*. DOI: 10.1101/2021.09.24.461498v3

---

For any further information, please contact: [Pablo Capilla-Lasheras](https://scholar.google.com/citations?hl=en&user=5JMTO-kAAAAJ&view_op=list_works&sortby=pubdate), email: pacapilla@gmail.com

## Code:

All R code is available in [`scripts`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/scripts). Scripts are numbered in order of use from 1 to 10. Scripts titles should be self-explanatory, but each script contains a description section with further information. The paths provided to import datasets into R assume your location is the main folder of the repository (i.e., the general folder `meta-analysis_variation_urban`). The R version used for this project was 4.2.0.

## Folders:

[`data`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/data): contains two files and one sub-folders: 
* [`META-ANALYSIS - URBAN_RURAL VARIATION - DATA BASE.xlsx`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/data): initial dataset with effect sizes collected from the literature.
* [`META-ANALYSIS - URBAN_RURAL VARIATION - DATA BASE - ReExtracted.xlsx`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/data): effect sizes from 10 studies, randomly chosen ([`List of study IDs data extracted from.csv`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/data/) contains the list of study IDs randomly generated), re-extracted to validate literature data extraction.
* [`processed_RDS_data_files`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/data/processed_RDS_data_files): processed dataset ready for analysis. Scripts in the folder [`scripts/1_data_processing`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/scripts/1_data_processing) generate these data tables.

[`scripts`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/scripts/1_data_processing): R code to produce results and figures included in the manuscript. It contains three subfolders:
* [`1_data_processing`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/scripts/1_data_processing): scripts to clean and generate tidy tables for further analysis.
* [`2_meta_models`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/scripts/2_meta-models): scripts for meta-analysis, meta-regressions and visualisation of results.
* [`R_library`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/scripts/R_library): help functions to run the rest of the scripts.

[`plots`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/plots): contains all figures (main and supplementary ones) created for this study, all of which were created using the R package [ggplot2 v.3.3.2](https://cran.r-project.org/web/packages/ggplot2/index.html).

[`models`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/models/): contains three subfolders with models included in [`Table S2`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/models/Table_S2) and [`Table S4`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/models/Table_S4) of the manuscript. A third sub-folder, [`land_cover_models`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/models/land_cover_models) contains all models whose results are presented in Figure 4 (i.e., meta-regressions including land cover variables as moderators).

[`literature_search`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/literature_search): contains the references found by our searches [`List of studies screened.xlsx`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/literature_search), which we performed in Web of Science and Scopus, after removing duplicates. The list contains details of which studies were inspected and which studies contained data that we extracted. This folder, file [`Benchmark papers.xlsx`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/main/literature_search), also contains the assessment of literature search comprehensiveness. Re-extracted effect sizes to check data extration can be found in [`Data_validatoin_reextraction_effect_sizes.xlsx`](https://github.com/PabloCapilla/meta-analysis_variation_urban/tree/master/literature_search).

## Notes

See details of the licence of this repository in [`LICENSE`](https://github.com/PabloCapilla/meta-analysis_variation_urban/blob/main/LICENSE).
