# MCF7 cell line reconstruction

## Running the case study

Please run all scripts with the human_ts_models repository folder as your Python root.

1. Run the *data_download* and *data_preprocessing* scripts inside the *omics* folder to download the CCLE dataset and 
preprocess the expression levels
2. Run the *reconstruction_with_task_evaluation* script to generate the models and evaluate metabolic tasks
3. Run the *task_heatmaps* script to generate task evaluation heatmaps.

## Missing steps

* Add simulation scripts
    * Sampling with hit-and-run or geometric approach
    * pFBA - might be unreliable