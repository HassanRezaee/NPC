---
title: "A Process Convolution Model for Crash Count Data on a Network"
author: "Hassan Rezaee"
date: "8/27/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Abstract

Crash data observed on a road network often exhibit spatial correlation due to unobserved effects with inherent spatial correlation following the structure of the road network. It is important to model this spatial correlation accounting for the road network structure. In this study we introduce the network process convolution (NPC) model. In this model, the spatial correlation among crash data is captured by a spatial random effect Gaussian Process (GP) constructed through a kernel convolution approach. The GP's covariance function is based on path distances computed between a limited set of knots and crash data points on the road network. The proposed model offers a straightforward approach to predict crash frequency at unobserved locations where covariates are available, and to interpolate the GP values anywhere on the network. Furthermore, the proposed model can be implemented in INLA, which offers a very efficient estimation procedure compared to MCMC sampling algorithms. We tested our model on synthetic data and a real crash data set from the city of Ottawa in Canada and we compared our approach to the proper Conditional Autoregressive (pCAR) and Poisson Regression (PR) models. The results of the study indicated that although pCAR has comparable fitting performance, our NPC models outperforms pCAR in a prediction mode at the unobserved locations. This can be seen from the lower mean absolute error rates between cross validated crash counts, latent variable values, fixed effect coefficients, as well as better interval scores particularly for singletons. We argue that NPC provides a more straightforward process to incorporate the network structure in road crash modeling and it offers improved predictive capability for crash data on a road network.


