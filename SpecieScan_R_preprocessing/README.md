# Pre-processing guide for SpecieScan in R

# Installation: 
* download R & RStudio @ https://posit.co/download/rstudio-desktop/

* download R pre-processing algorithm for SpecieScan @ https://github.com/mesve/SpecieScan

# License: MIT

# Credits
* Cite SpecieScan algorithm: DOI: 10.5281/zenodo.8055426 and paper Végh & Douka 2023. SpecieScan: semi-automated taxonomic identification of bone collagen peptides from MALDI-ToF-MS. Bioinformatics 40(2), 1-12.


* The Preprocessing code was adapted and subsequently changed from the following code and paper: Maria Codlin, & Richter, K. K. (2022). ZooMS processing and clustering workflow using MALDIquant. Zenodo
Codlin MC, Douka K, Richter KK (2022). An application of zooms to identify archaeological avian fauna from Teotihuacan, Mexico. Journal of Archaeological Science;148:105692.
![image](https://github.com/mesve/SpecieScan/assets/25906714/12ccc9fc-b1f8-4b2f-a34b-bc69fec8b695)


# Usage
* Specific parameter choices may depend on the characteristics of your data. It is essential to consider your data's unique properties. The common parameters in the code are explained below:

## SNR (signal-to-noise ratio)

SNR is a measure of the signal's strength compared to noise level, used to distinguish peaks from background noise. The default SNR parameter used in mMass is 3.0 (hence the default in the R code), but this needs to be adjusted accordingly after observing the spectra using the following code for visualisation. Higher values will filter out more noise but may miss weaker peaks. If a lower value is chosen, it might detect more peaks but also more noise.Thus, this parameter will depend on your spectra.


```{r}

noise <- estimateNoise(spectra[[1]])
plot(spectra[[1]], xlim=c(1150, 1250))
ylim(c(0, 0.01)) # set y-axis limits to 0-0.01
lines(noise, col="red")
lines(noise[,1], noise[, 2]*3, col="blue") # here the *3 shows SNR of 3
for (i in 1:10){
  plot(spectra[[i]], xlim=c(3090, 3100))
  ylim(c(0, 0.01)) # set y-axis limits to 0-0.01
  lines(noise, col="red")
  lines(noise[,1], noise[, 2]*5, col="blue") # here the *5 shows SNR of 5
}

```

## Half-Window Size

Smoothing is used to reduce noise in the spectra and helps eliminate the baseline signal that can obscure peptide peaks. This parameter determined the width of the smoothing window. In the default code included, a half-window size of 100 was chosen with the SNIP method because it captured broader features in the data while reducing high-frequency noise. Adjust this parameter accordingly, which will affect the degree of smoothing. In our data it was found that a half-window size of 100 was a good balance between preserving underlying features of the data and reducing noise throughout our experimentation. It should be noted that this parameter will also depend on the instrument used. A larger half-window size will cover a wider range of data points when applying the smoothing algorithm. It is more effective at noise reduction because it averages data points over a larger span. If the spectra is very noisy, a larger half-window size is advised because it will reduce high-frequency noise that can obscure genuine peaks in the spectra. However, this way there is a larger chance of over-smoothing the data resulting in the loss of important spectral features (i.e., small but significant peaks might become overly blur). 

### Sensitivity Analysis

* Variability Assessment: We began by assessing the variability of baseline artifacts across our MALDI-ToF-MS spectra. We observed that baseline characteristics, including intensity and shape, varied across different samples.

* Iterative Testing: To determine the optimal number of iterations for baseline removal, we performed baseline removal with different iteration counts, ranging from 10 to 200. For each iteration count, we closely examined the impact on baseline removal quality and spectral integrity.

* Quality Metrics: To quantitatively assess the quality of baseline removal, we used two key metrics:
- **Signal-to-Noise Ratio (SNR)**: A higher SNR indicates better removal of baseline noise.
- **Root Mean Square Error (RMSE)**: A lower RMSE signifies improved baseline subtraction with less distortion of true signal peaks.

#### Parameter Tuning Explanation

- **Overfitting and Underfitting**: We carefully considered the trade-off between overfitting and underfitting during the baseline removal process. Using too few iterations can result in underfitting, leaving residual baseline artifacts. Conversely, using too many iterations may lead to overfitting, causing distortion of the true signal in our spectra.

- **Balance and Stability**: After conducting testing, we found that using 100 iterations struck a balance between effective baseline removal and maintaining spectral stability. Beyond this point, increasing the number of iterations did not significantly improve baseline removal quality but introduced instability and over-smoothing in the spectra.

- **Computational Efficiency**: Practical considerations also played a role in our parameter choice. While more iterations might theoretically provide better baseline removal, it can be computationally expensive. Using 100 iterations allowed us to achieve our desired baseline removal quality without imposing excessive computational overhead.

- **Reproducibility**: We ensured that our choice of 100 iterations yielded consistent and reproducible results across various samples and datasets.

* Based on our sensitivity analysis and parameter tuning, we chose to use 100 iterations for the baseline removal step in our MALDI-ToF-MS data preprocessing pipeline. This parameter choice provides effective baseline removal while maintaining spectral integrity, computational efficiency, and reproducibility across our datasets.

Section of code in question:

```{r}
spectra <- removeBaseline(spectra, method = "SNIP", iterations = 100)
}

```
## Monoisotopic Peaks Parameters
It allows for the identification and quantification of peptide fragments accurately. The `monoisotopicPeaks` function plays a central role in this process. Here, we explain our choices of parameters for this function and provide insights from our parameter testing and sensitivity analysis.

### Parameter Settings

In our analysis pipeline, we used the following parameters for monoisotopic peak detection:

- **minCor**: 0.95
- **tolerance**: 1e-4
- **distance**: 1.00235
- **size**: 2L:10L

We tested these and other parameters using sensitivity analysis, during which we sought to optimise the detection of monoisotopic peaks while minimising false positives and false negatives. 

Other parameter values tested and have the potential to be used in different types of spectra/instrument:
- **minCor**: 0.90, 0.92, 0.94 *Lowering the 'minCor' threshold allows for more relaxed correlation requirements between peaks within an isotopic cluster. This is useful when you have noisy data or when you want to capture weaker peaks.

- **tolerance**: 1e-3, 5e-5, 1e-5 *Adjusting this value can change the allowed mass difference between peaks in an isotopic cluster. A lower 'tolerance' value enforces stricter mass matching criteria.

- **distance**: 1.003, 1.001, 1.004 *Modifying this value (isotopic distance) changes the expected mass difference between isotopic peaks. A higher 'distance' value allows for a wider range of isotopic distances between peaks.

- **size**: 3L:15L, 4L:12L, 2L:8L *Adjusting this parameter changes the allowable range of cluster sizes. A wider range (e.g., 3L:15L) accommodates larger isotopic clusters, while a narrower range (e.g., 2L:8L) focuses on smaller clusters. In Végh & Douka 2023. the specific choice was size=2L:10L, which accommodated potential factors, like partial isotope distributions of overlapping peptides and the impact of deamidation on peak patterns. Considering the potential impact of deamidation, which can lead to variations in peak locations and patterns, this setting strikes a balance between sensitivity and specificity, enabling the detection of peaks influenced by deamidation while maining a level of moderate stringency to minimise false positives. 

#### Sensitivity Analysis

- We started with default parameter settings for the `monoisotopicPeaks` function but found that they required further refinement for our spectra.

- A sensitivity analysis was conducted by systematically varying `minCor`, `tolerance`, `distance`, and `size`. I evaluated the impact of these parameter variations on the number of true positive, false positive, and false negative detections.

#### Optimization Criteria

For optimisation we aimed to strike a balance between the accurate detection of monoisotopic peaks (high recall) and minimising incorrect detections (high precision).

#### Parameter Choices

After testing, we arrived at the given parameters, but the other listed parameters can also be used with different types of spectra. These settings consistently provided a good trade-off between sensitivity and precision in detecting monoisotopic peaks across a range of spectra in our dataset.We encourage users to consider their specific data characteristics and analysis goals when choosing parameters for monoisotopic peak detection. The settings we've provided here represent a starting point and can be adapted to suit individual dataset requirements.


## Peak Picking Method (SuperSmoother)

Different peak picking algorithms have different strengths and weaknesses. 'Supersmoother' was chosen based on its performance at retaining peptide peaks and suppressing noise, which influences the number and quality of detected peaks. 

- **Method**: SuperSmoother

### Parameter Settings

For peak picking using the SuperSmoother method, we used the following parameter settings:

- **Half-Window Size**: 20. Through sensitivity analysis and visual inspection, we found that this half-window size effectively reduced high-frequency noise while maintaining the integrity of important spectral features. Smaller window sizes resulted in less noise reduction, while larger sizes over-smoothed the data.

- **Signal-to-Noise Ratio (SNR)**: changes based on the appearance of the peaks and noise. A higher SNR threshold reduces the likelihood of false positive detections.



