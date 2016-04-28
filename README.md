This repository contains the data and code necessary to reproduce the analysis in the manuscript "Bat activity during autumn relates to atmospheric conditions: implications for wind energy development" in revision at the *Journal of Mammalogy*.  

Authors (affiliation):
- Adam D. Smith<sup>*</sup> (Department of Natural Resources Science, University of Rhode Island, 105 Coastal Institute, 1 Greenhouse Road, Kingston, RI 02881 USA)
- Scott R. McWilliams (Department of Natural Resources Science, University of Rhode Island, 105 Coastal Institute, 1 Greenhouse Road, Kingston, RI 02881 USA )

<sup>*</sup>corresponding author; current address: U.S. Fish & Wildlife Service, National Wildlife Refuge System, Southeast Inventory & Monitoring Branch, 135 Phoenix Road, Athens, GA 30605 USA)

The easiest way to access the data and code necessary to reproduce the entire analysis is to [fork and clone](https://help.github.com/articles/fork-a-repo/) this repository.

To see an overview of the analysis (Supporting Information S2 of the manuscript):

1. Fork this repository
  - all code and data is available in this repository
2. Open [R](http://www.r-project.org) on a computer with internet access
3. Install the [devtools](http://cran.r-project.org/package=devtools) and [rmarkdown](http://cran.r-project.org/package=rmarkdown) packages, if necessary 
  - install.packages("devtools")
  - install.packages("rmarkdown")
4. Render the file to a word document
  - rmarkdown::render('PATH/TO/Bat_migration_acoustics.Rmd', output_format = 'word_document')
5. Open the resulting Word document
  - we've included this file (Bat_migration_acoustics.docx) if you'd rather download it and open it directly

