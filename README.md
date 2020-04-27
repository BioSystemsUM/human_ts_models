# human_ts_models
A repository for tissue-specific model reconstructions using troppo. 
Intended to separate package code with model creation and data analysis scripts. 
Ask to become a collaborator if you want to contribute.

## Project structure
There are two main folders in this project:
* **projects** - Hosts individual projects for each case study
* **shared** - Contains resources to be made available for all projects

## Creating your own project
Ask first to create a branch for yourself and work on that branch. This helps 
Create a folder inside the *projects* folder of this project with your desired 
model/project name. Your data (if it's small enough to fit GitHub) should be 
inside this folder. A recommendation for this is to organise your folder with 
the following structure:
* **data/** - Put your omics data here or any other files related to your specific 
case study. Do not commit large files! This folder can actually be empty, as long 
as your scripts are able to download the data from the original sources.
* **scripts/** - Python (or otherwise) scripts that can be run individually for 
anything related with your analysis. IPython notebooks should also be stored in 
this folder (or in a **notebooks/** folder, if you prefer)
* **results/** - File outputs from the developed scripts 
(results from individual runs of your scripts). It is recommended that you create 
sub-folders if you're running scripts for multiple conditions/datasets or if you're 
trying various parameters for your algorithms.
* **support/** - Any file that isn't data from your case study but is needed for 
your specific project. This includes intermediate files you might have generated 
(such as a simplified metabolic model tailored to suit your needs or a pickled 
object to speed up your scripts)
* **src/** - Source-code that isn't dependent on anything from your own project 
but can be used by your scripts.

Please do not commit any changes from someone else's project unless 
you've been told to do so.

Once you feel your project is matured enough to share with everyone, you can
submit a pull request and make your changes available in the master branch

## Contributing with shared resources
If you find that a shared resource could be complemented with your own code/modifications,
submit a pull request and let other people review your changes first.