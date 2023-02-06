# Supervised Classes, Unsupervised Mixing Proportions: Detection of Bots in a Likert-Type Questionnaire

Michael John Ilagan and Carl F. Falk

This repository contains computer code for the article with the following citation:  
Ilagan, M.J. & Falk, C.F. (2022). Supervised classes, unsupervised mixing proportions: detection of bots in a Likert-type questionnaire. Educational and Psychological Measurement. [https://doi.org/10.1177/00131644221104220](https://doi.org/10.1177/00131644221104220)

## Abstract

Administering Likert-type questionnaires to online samples risks contamination of the data by malicious computer-generated random responses, also known as bots. 
Although nonresponsivity indices (NRIs) such as person-total correlations or Mahalanobis distance have shown great promise to detect bots, universal cutoff values are elusive. 
An initial calibration sample constructed via stratified sampling of bots and humans---real or simulated under a measurement model---has been used to empirically choose cutoffs with a high nominal specificity. 
However, a high-specificity cutoff is less accurate when the target sample has a high contamination rate.
In the present article, we propose the supervised classes, unsupervised mixing proportions (SCUMP) algorithm that chooses a cutoff to maximize accuracy. 
SCUMP uses a Gaussian mixture model to estimate, unsupervised, the contamination rate in the sample of interest. 
A simulation study found that, in the absence of model misspecification on the bots, our cutoffs maintained accuracy across varying contamination rates.

## Contents

This repository contains the following files.

* `helper-general.R`, an R script defining functions for computing nonresponsivity indices, plotting, and evaluating classification accuracy.
* `helper-superv.R`, an R script defining functions implementing SCUMP.
* `sim.R`, an R script carrying out the simulation study reported in the article.
* `figs.R`, an R script for producing the Figures in the article.

## Data

The raw data used in the simulation study is answers to the Humor Styles Questionnaire from the Open Psychometrics Project.  
[http://openpsychometrics.org/_rawdata/HSQ.zip](http://openpsychometrics.org/_rawdata/HSQ.zip)

## Instructions

To reproduce the results in the article, have all the R scripts in your working directory as well as the data file `data.csv` from the Open Psychometrics Project dataset.
Run `sim.R` to produce the R workspace of the results.
Then run `figs.R` to produce the Figures in PDF format.
