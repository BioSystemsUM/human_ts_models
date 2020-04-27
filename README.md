# human_ts_models
A repository for tissue-specific model reconstructions using troppo. Intended to separate package code with model creation and data analysis scripts. Ask to become a collaborator if you want to contribute.

## Contributing to this project

Create a folder in the root of this project with your desired model/project name. Your data (if it's small enough to fit GitHub) should be inside this folder. A recommendation for this is to organise your folder with the following structure:
* **data/** - Put your omics data here or any other files related to your specific case study. Do not commit large files!
* **scripts/** - Python (or otherwise) scripts that can be run individually for anything related with your analysis. IPython notebooks should also be stored in this folder (or in a **notebooks/** folder, if you prefer)
* **results/** - File outputs from the developed scripts (results from individual runs of your scripts). It is recommended that you create sub-folders if you're running scripts for multiple conditions/datasets or if you're trying various parameters for your algorithms.
* **support/** - Any file that isn't data from your case study but is needed for your specific project. This includes intermediate files you might have generated (such as a simplified metabolic model tailored to suit your needs or a pickled object to speed up your scripts)
* **src/** - Source-code that isn't dependent on anything from your own project but can be used by your scripts.

Please do not commit any changes from someone else's project unless you've been told to do so.
