# MCF7 cell line reconstruction

## Running the case study

Please run all scripts with the human_ts_models repository folder as your Python root.

1. Run the *data_download* and *data_preprocessing* scripts inside the *omics* folder to download the CCLE dataset and 
preprocess the expression levels
2. Run the *reconstruction_pipeline_comparison* script to generate models for various preprocessing threshold
combinations
3. Run the *task_heatmaps* script to generate task evaluation heatmaps.

## Missing steps

* Add simulation scripts
    * Sampling with hit-and-run or geometric approach?
    * pFBA - might be unreliable
    * FVA
    
* Reconstruct normal breast tissue with GTEx