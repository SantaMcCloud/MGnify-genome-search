# MGnify genome search

This tool is a copy of from the [MGnify Data analysis notebooks](https://shiny-portal.embl.de/shinyapps/app/06_mgnify-notebook-lab?jlpath=mgnify-examples/home.ipynb). 

The goal is to upload fasta files which then will get queued against the MGnify database if they exit there or not and if which taxa are they. For this for each MAG a sourmash sketch will be created and will be used for the whole process to find the novelty and the classification via MGnify.

## Inputs

The tool only has one input, which is the path to the directory where all the fasta files are stored. 

## Outputs

There are 3 outputs which are the followed:

- 1 table where the matches of the MGnify catalogues are stored for each MAG which got a match. Together with this there are more information stored!
- 1 table where each MAG are classified by the best match from MGnify.
- 1 table where the novelty is stored. Therefore if the MAG is novel the result in the second column is true otherwise it is false.
