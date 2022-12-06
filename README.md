# GenoMine
Shiny app for data exploration of RNAseq data. (work in progress)

## Home:

The application is splited in 4 main steps:

* First: upload your expression data and clinical info data and explore your data.
* Second: Select your normalitzation model and explore clustering using multiple variables from your clinical info.
* Third: Perform a diferential expression analysis using diferent models.
* Fourth: Pathway analysis using diferent GO from yout DE results.

![My Image](images/menu.png)

## Data:

In this section the user can upload their experssion data in raw counts and the clinical or variable data from their samples.
For the counts expression files, genes must be on the rows and samples on columns.
Clinical/variable info data must be samples on the rows and variables on columns.
Please select if your data have headers or not on the selector button.

After uploading the data a summery of your data will appear on Summary Info section, a table with your sample info can be visualize in the Table info secction and the table with your counts can be visualize in Counts section.

![My Image](images/data.png)

![My Image](images/clustering.png)

![My Image](images/DE.png)

![My Image](images/pathway.png)
