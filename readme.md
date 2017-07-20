
# C-PAC addons

Here are some notes and tools relating to using [C-PAC](https://github.com/FCP-INDI/C-PAC).

* `make_sca_report.py` - A little tool that takes a seed-based correlation analysis (SCA) output folder and produces an overview of the results in a self-contained HTML file. Usage: `python make_sca_report.py <CPAC_gpa_config.yml>` (where `CPAC_gpa_config.yaml` is the CPAC configuration file for the group analysis), which will parse the corresponding output directory and will create an HTML file with the results in the parent directory of the output directory.






# Some notes about using C-PAC

## Building a subject list

```
cpac_setup.py /path/to/data_config
```


## Setting up contrasts
You can make your own contrasts file, and F-tests can be added by adding columns, but it is important to note that they should be called `f_test_<SOMETHING>`, otherwise they are interpreted as columns from the design matrix (which do not exist).


## Running

To run, use:

```
$ bash                   # will ensure that you use bash
$ source activate cpac   # activates the CPAC conda environment
```

You can build a subject list using:
```
$ cpac_setup.py data_config.yml
```

You can then run an analysis using:
```
$ cpac_run.py /path/to/pipeline_config.yml /path/to/CPAC_subject_list.yml
```

However, I believe this runs a subject-level analysis only.
For the group level:
``` 
$ ipython
# import CPAC
# CPAC.pipeline.cpac_group_runner.run(config_file,pipeline_output_folder)
```

e.g.
```
#  CPAC.pipeline.cpac_group_runner.run('/brains/pianists/pipeline_config_pianists-3mm-fsl.yml','/brains/pianists/output/pipeline_pianists-3mm-fsl')
```





## Understanding CPAC group level output structure for SCA

In your output folder, you have nested folders corresponding to each analysis step, e.g.
`group_analysis_<<NAME>>/group_model_<<NAME>>/sca_roi_files_to_standard_smooth_fisher_zstd/REST/_compcor_ncomponents_5_selector_pc10.linear0.wm1.global0.motion1.quadratic0.gm0.compcor1.csf1/_bandpass_freqs_0.01.0.1/_mask_rois_1mm/_fwhm_6`

Corresponds to the seed-based correlation analysis (`SCA`) applied to the scan called (`REST`), applying various nuisance regressors (`compcor_..._wm1.global...`) then bandpass filtering (`bandpass_freqs_0.01.0.1`) and then applying the ROI mask (`_mask_rois_1mm`) and finally smoothing with a 6mm kernel (`_fwhm_6`).

In this folder, you will see one subfolder for every seed. If you had defined four seeds, you will see four folders like this:
```
_fisher_z_score0
_fisher_z_score1
_fisher_z_score2
_fisher_z_score3
```

Inside each folder, you will see a subfolder called `sca_ROI_X` where `X` is the number of the seed (again). That brings you to the main results folder for that seed, which contains four subdirectories:

* `merged` - contains an output image with one volume for each subject, giving the connectivity map for that subject and the given seed.
* `model_files` - files pertaining to the model, i.e. contrasts, f-tests, design matrices, voxel time servies, etc. 
* `rendered` - rendered group z-maps for each of the contrasts and f-tests you performed, cluster-corrected. by convention, maps that are empty are drawn with a yellow haze over it, so that you don't need to look too closely to see that you're out of luck ;)
* `stats/clusterMap` - text file list of all clusters, their size, extent, p-value, etc.
* `stats/threshold` - thresholded z-maps
* `stats/unthreshold` - unthresholded z-maps





