## Miscellaneous scripts, processing pipelines and analysis steps for colleagues

#### annotate-maf-wrapper.R
example script, where annotate-maf-wrapper.R is located in `${HOME}`   
* `Rscript <script_location> -m <maf_location> -f <facets_output.Rdata_location> -o <output_directory>`   
* `Rscript /juno/home/kreitzec/annotate-maf-wrapper.R -m /juno/work/schultz/subhi/mafs_for_mafAnno/F12-6-p1--F12-blood.union.v4.annotated.maf -f /juno/work/subhi/Ingo_WES_Facets/s_F12-6-p1.gz.Rdata -o /juno/work/schultz/subhi/F12-6-p1_annotated.txt`


## Data access various formats and important links (MSK related); VPN required    
[MSK-knowledgeSystems](https://github.mskcc.org/)  access only via VPN   
[Delphi-Sample-Tracker](https://delphi.mskcc.org/sample-tracker/home)   
[High Performance Computing](http://mskcchpc.org/) Questions related to juno and server access   
[CVR Molecular Diagnostic service](https://cvr.mskcc.org:8083/search) - open @Chrome: mskcc Christoph    
[DMP-2022 DMP Data Directory](https://github.mskcc.org/knowledgesystems/dmp-2022) #' Page 404 not found / need to login @Enterprise (msk initials)    
[OncoKB annotated Mutations](https://github.mskcc.org/knowledgesystems/oncokb-annotated-msk-impact) daily updates with OncoKB annotations / need login as well    
[Darwin](https://ddp.mskcc.org/search) Clinical Interpretation tool; all clinical data associated with tumor specimen; specifically pathology and radiology events and treatment data.   


## Survival analysis; be cautions about left-truncation and length bias
Read this handy [article](https://towardsdatascience.com/how-well-do-you-really-know-survival-analysis-221885b89a8e#:~:text=Left%20truncation%20occurs%20when%20data,study)    
Use this [Github](https://github.com/slb2240/delayed_entry_clin_genom_studies/blob/main/crc_stage_iv_os_dx.R) for R   
And some other useful [publication](https://pubmed.ncbi.nlm.nih.gov/34734967/)   
