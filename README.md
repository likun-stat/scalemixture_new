# SCALE MIXTURE for Flexible Extremal Dependence Model

## Research ongoing with Dr. Ben Shaby

[slides from JSM] (https://drive.google.com/file/d/1xUYZZ9jN4xxxkg57UQWvTDuV3Jp6d4p2/view?usp=sharing) .

We try to model dependence in extremes of spatial processes: 

![equation](http://latex.codecogs.com/gif.latex?%5C%7BX%28s%29%3Bs%5Cin%5Cmathcal%7BS%7D%5Csubset%5Cmathbb%7BR%7D%5E2%5C%7D)

### 1.  Flexible Spatial Modeling
We wish to fit spatial models encompassing both AD and AI, and let the data choose in between.

![equation](http://latex.codecogs.com/gif.latex?X%28s%29%3DRW%28s%29&plus;%5Cepsilon%28s%29%2C%20%5Cepsilon%28s%29%5Cstackrel%7Biid%7D%7B%5Csim%7D%20N%280%2C%5Ctau%5E2%29)

where

![equation](http://latex.codecogs.com/gif.latex?R%7C%5Cdelta%5Csim%20Pareto%28%5Cfrac%7B1-%5Cdelta%7D%7B%5Cdelta%7D%29)

and 

![equation](http://latex.codecogs.com/gif.latex?W%28s%29%3DF%5E%7B-1%7D%28%5CPhi%28Z%28s%29%29%29%3DT%28Z%28s%29%29),
![equation](http://latex.codecogs.com/gif.latex?Z%28s%29%5Csim%20N%280%2C%5CSigma%28%5Clambda%2C%5Cgamma%29%29).

### 2. Inference for Scale Mixtures

We are restricting focus to the tail by censoring, so that the low values won't affect the inference of the extremal dependence structure:

![plot1](www/6.png)

We sample the latent process as follows:

- Draw the smooth process _X*(s)_ given the noisy process _X(s)_.

(No truncation!)

- Draw the noise given the smooth process _X*(s)_.

(Truncated, but univariate!)

![gif](www/anime.gif)
