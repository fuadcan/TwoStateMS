
# Two State MS-ARFIMA 

This repository consists R codes and helper files for Two State Markov Switching ARFIMA (2MS-ARFIMA) estimation proposed by Tsay and Härdle (2009). The original codes by Tsay and Hardle (2009) were translated from GAUSS to R and some errors are corrected. The repo also includes output gap analysis by [Beylunioglu et.al. (2018)](https://www.degruyter.com/view/journals/snde/22/3/article-20170043.xml).


## Outline
The main code fits the model is `MSArfimaFit.R` with two inputs; series to be fit and type of the fit. The type can be either DMS,DM,DS or D that controls the parameters to state switch: for DMS, $(d,\mu,\sigma)$ switche between states whereas for D only $d$ switches so on.

`MSArfimaFit.R` calls `lnviDMS`,`lnviDM`,`lnviDS` or `lnviD` with respect to the type given. These codes calculates likelihoods for a given series and a given parameter vector ($d_1,d_2,P_{11},P_{22},\mu_1,\mu_2, \sigma_1, \sigma_2$) where in $\mu_1 = \mu_2, \sigma_1 = \sigma_2$ are assumed in `lnviD2.R`. The parameter vector is given to `lnviD*` exogenously so that `MSArfima.fit` optimizes the function over the vector. The initial parameters and boundary values for optimization is given inside MSArfimaFit.R file and three different cases of initial vectors are written for each type, if the optimization does not converge, the rest is referred. The initial sets are chosen to meet the requirements for output gap analysis in Beylunioglu et.al (2018). Be sure to adjust them your analysis.

The rest of the codes support the analysis of Beylunioglu et.al (2018) which can be decomposed to three steps; estimating parameters $d,P,\mu,\sigma$ etc., reporting descriptive statistics of the results and generating plots. All outputs are saved to disk.

Codes are designed in R project version 3.41 for Windows and developed in both Ubuntu 16.04 and Windows 8 64bit. To run the code trouble-free, please check your R version and save clone the repository into Documents folder.

## Usage
Run all the commands in main.R. The file contains commands used for reading, estimating parameters, reporting and plotting results for all pairs for each data. The preliminary results are saved in Results folder as .rda files. NOTE: The tables in these files are re-formed for the paper but both shares same outline.

### Codes
- lnviDMS2.R, lnviDM2.R, lnviDS2.R, lnviD2.R: Estimates state switching parameters (see above) and transition probabilities for a given series.

    Inputs:
    1. b: Initial set of parameters
    2. w: series to be fitted MS-ARFIMA$ (0,d,0) 
		 
    Output: Log-likelihood of the estimate


	 
- convDLV.R: Calls lnviD**2.R for all log-differentiated gap series.
	
    Input:
		
			
	1. yearOrRegion: Name of the data (1930,1940,Europe+G7,Europe+S&P,G7+S&P)
		
	Outputs: Parameter estimates for all gap series (both saved to disk and returned)
		


- plotAll.R: Calls results which are corrected in main.R, extracts the path by using dlvPath.R, plots the path and saves to pathplots folder.
	
	Input:
		
	1. yearOrRegion: Name of the data (1930,1940,Europe+G7,Europe+S\&P,G7+S\&P)
		
	Output: Plots are saved to disk
	


- plotRejs.R: Calls results which are corrected in main.R, classifies results with $ d<1 $ as S (for stationary or mean reverting) and $ d\geq 1 $ as N (for non-stationary), calculates and plots the ratio of the pairs with $ d<1 $ by time for a given dataset; then saves to plots folder.
		
    Input:
			
    1. yearOrRegion: Name of the data (1930,1940,Europe+G7,Europe+S&P,G7+S&P)
			
    Output: Plots are saved to disk	


- main.R: Calls all above codes. Additionally processes outputs and saves to disk. The table in paper is generated and saved into ReplicationFiles folder as table.csv.


## References

Tsay, Wen-Jen, and Wolfgang Karl Härdle. "A generalized ARFIMA process with Markov-switching fractional differencing parameter." Journal of Statistical Computation and Simulation 79.5 (2009): 731-745.

Beylunioglu, Stengos, Yazgan. "Regime Switching Output Convergence" SNDE (2018), In publishing
