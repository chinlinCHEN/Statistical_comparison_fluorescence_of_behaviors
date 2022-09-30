# Statistical_comparison_fluorescence_of_behaviors

This repository is part of the analysis pipeline of https://github.com/NeLy-EPFL/Ascending_neuron_screen_analysis_pipeline/ for the manuscripts [**Ascending neurons convey behavioral state to integrative sensory and action selection centers in the brain**] in bioRxiv (https://www.biorxiv.org/content/10.1101/2022.02.09.479566v1), which focuses on comparing the fluroscent neural actiivty during different behavior epochs.

First, a statistical analysis on the data are retrieved from a distributed data structure instead of a centralized dataframe to maintain the user-friendly daily check on each individual target's results. 

Second, because the nature of large-sampling time-series data causes false positive result, this analysis uses bootstrapping and voting strategies to avoid oversized sample numbers and thus keep the efficacy of statistical hypothesis to comparison. 

Note: this repository is to archived the code for easier understanding. It won't work with the proper data structure due to the large size of the raw data. The intact data structure should refer to https://github.com/NeLy-EPFL/Ascending_neuron_screen_analysis_pipeline/. 

## Content
- [Installation](#installation)
- [Steps](#steps)
  - [Overview](#overview)
  - [Statistical comparison with bootstrapping and voting approahces](#statistical-comparison-with-bootstrapping-and-voting-approahces)
  - [Result visualization](#result-visualization)
  
## Installation
To install the AN environment for running Python scripts, please refer to https://github.com/NeLy-EPFL/Ascending_neuron_screen_analysis_pipeline/

## Steps

### Overview
#### Computation end
```1-statistics_dFF_comparison_whole dataset.py``` does:

(a). It retrieves time-series data set (```.dic```) of each trials of experiments across folders of 50 genotypes.

(b). It overlays the neural activity signals during the same type of behacvioral epochs and iterate this action through differnet types of behavior.

(c). It perform ANOVA and Tukey posthoc comparison and output the datafram summarize statistical result (```.pkl```) and the overlaid traces of each behavior period.

```2-plot_matrix.py```

(a). It reads the datafram summarize statistical result (```.pkl```) and output a matrix plot.

#### User end
one can navigate the data structure with graphical user interface to check the intermediate or final analysis results together with other input such as images and neural morphology.

<p align="center">
  <img align="center" width="780" src="/images/statistic_analysis-data_retrieval_diagram-01.png">
</p>

### Statistical comparison with bootstrapping and voting approahces
1. A baseline of neural activity trace is derived by thresholding with otsu filter if the Shapiro-Wilk normality test of the datapoint distribution reject the null-hypthesis, which indicates that the trace is not background noise but with signals. 

2. ANOVA and Tukey comparison are performed on each group of bootstrapped samples from the datapoint of baseline and different behavioral epochs. This process is iterated 6 times and the significant difference against baseline are used for voting to determine if the cell significantly activate during each behavioral epochs. The significance test would pass if the voting is higher than 50%. 

3. In the meantime, this script also output the mean activity and 95%CI during different behavioral epochs and the pickle file of statistical summary (```.pkl```).

<p align="middle">
  <img align="middle" width="400" src="/images/statistics_analysis_approach-01.png">
  <img align="middle" width="600" src="/images/20190311_SS27485-fly2-posthoc_compar_dFFofBeh_ROI_2.png">
</p>

### Result visualization

Finally, ```2-plot_matrix.py``` reads the pickle file of statistical summary (```.pkl```) to output he matrix plots.
<p align="center">
  <img align="center" width="900" src="/images/dFF_comparison_statistics_whole_dataset.png">
</p>








