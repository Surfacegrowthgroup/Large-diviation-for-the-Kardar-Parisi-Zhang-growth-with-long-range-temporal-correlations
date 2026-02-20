# Simulation codes and data for "High-precision distributions for kinetic roughening in the Kardar-Parisi-Zhang growth with long-range temporal correlations"

## Contents

There are three parts matlab code and data of the paper. Simulations and results of the paper can be obtained by these files. Every code file has an introduction in the beginning of the code, and some information is as fllows.

### Running environment


All of three codes work with Matlab R2024b and the Add-ons including Parallel Computing Toolbox and Statistics and Machine Learning Toolbox. 



The five code files represent the large deviation algorithm for KPZ and EW growth, and the error bars for the skewness and kurtosis, respectively. Some critical parameters and their meaning are as follows.

| Parameter        | Meaning                                |
| ---------------- | -------------------------------------- |
| dt               | the time difference interval           |
| nabla            | the nonlinear parameter                |
| L                | the system size                        |
| T                | the iteration times                    |
| thc              | the temporally correlated exponent     |
| thd              | the biased parameter                   |
| P                | the replacement rate of random numbers |
| dA               | the acceptable deviation               |

### The large deviation part of KPZ or EW

In order to obtain the high-precision height distribution of the KPZ and EW equation subject to temporal correlations, the process of these code is as follows.
 1.Determining the biased parameter
 2.Calculating the distribution under the chosen parameter
 3.Re-weighting the histogram and renormalization

### The error bar part of KPZ or EW

The code can obtain the height samples using simple sampling. If we want the error bars of skewness or kurtosis, we can calculate some samples of particular height data, then calculate the standard deviation.
