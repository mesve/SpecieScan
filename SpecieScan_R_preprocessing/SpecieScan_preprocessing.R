# Preprocessing for SpecieScan 

## SpecieScan: Automated Species Identification Algorithm for Bone Fragments from MALDI-ToF MS Spectra 

##v0.2.0

## For usage and guidance see README.md and comments within this script. This preprocessing script is to preprocess MALDI-ToF MS data for automated species identification.

## License: MIT

## Cite paper: Végh & Douka 2023. SpecieScan: Automated Species Identification Algorithm for Bone Fragments from MALDI-ToF MS. Bioinformatics

## & cite Github DOI: 10.5281/zenodo.8055426

## Authors: Emese Végh

## Contact: Emese Végh[mesve3@gmail.com]

## The following Data Preprocessing code was adapted from 'Codlin, M. C., & Richter, K. K. (2022). Data for ZooMS analysis of avian fauna from Teotihuacan, Mexico, for Codlin et al. 2022 [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6363114
Codlin MC, Douka K, Richter KK (2022)' and 'An application of zooms to identify archaeological avian fauna from Teotihuacan, Mexico. Journal of Archaeological Science;148:105692.'
and further developed here. 

```

## Load packages
```{r, include=FALSE}
library(factoextra)
library("MALDIquant")
library("MALDIquantForeign")

```


# IMPORT SPECTRA

## Import spectra straight from the MALDI:
### Windows:
spectra <- importBrukerFlex("D:\\Documents\\MALDI_FOLDER", verbose = FALSE)
### Mac:
spectra <- importBrukerFlex("~/Documents/MALDI_FOLDER", verbose = FALSE) #Delete "User/user_name" before the folder it was saved in 

# EXPORT mzML to rename files
## Windows:
exportMzMl(spectra, path = "D:\\Documents\\MALDI_FOLDER-Converted")
## Mac:
exportMzMl(spectra, path = "~~/Documents/MALDI_FOLDER-Converted")

# RENAME TECHNICAL REPLICATES IN YOUR FOLDER SO THEY ALL HAVE A UNIQUE IDENTIFIER (EITHER MANUALLY ONE-BY-ONE OR USE PYTHON CODE IN GIT REPO "RENAME_FILES") -- IMPORT THE MZML OR TXT DATA WITH A UNIQUE IDENTIFIER FOR EACH SAMPLE


# IMPORT RENAMED FILES

```{r,include=FALSE}
## Windows:

spectra<-importMzMl("D:\\Documents\\MALDI_FOLDER-Converted")

## Mac:
spectra<-importMzMl("~/Documents/MALDI-FOLDER-Converted")

```

# Save workspace to load later -- optional 

```{r}
save.image(file="DATA_r.RData")
```

#SPECTRA PROCESSING

##Make plots to visualise spectra:
```{r}

for (i in 1:10){ #to visualise the 1st to 10th spectra -- change 10 to a higher number if you want to see more 
  plot(spectra[[i]])
  plot(spectra[[i]], xlim=c(1500, 1600)) #xlim sets the limits of the x-axis. Edit to zoom in on regions of interest
}
```
##Quality control

```{r}
any(sapply(spectra, isEmpty))
table(sapply(spectra, length))
all(sapply(spectra, isRegular))
```
##Pre-processing
```{r}
spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)

spectra <- removeBaseline(spectra, method="SNIP", iterations=100)

spectra <- calibrateIntensity(spectra, method="TIC")
```

## Determine noise to apply SNR ratio

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
#PEAK PICKING OF ALL TECHNICAL REPLICATES
```{r}
peaks <- detectPeaks(spectra, method="SuperSmoother", halfWindowSize=20, SNR=3.0) # Change SNR to a number that excluded the noise on the plots
peaks <- monoisotopicPeaks(peaks, minCor=0.95, tolerance=1e-4, distance=1.00235, size=2L:10L) # Adjust parameter values -- see README.md in the SpecieScan R repo
```
## View picked and deisotoped peaks

```{r}
# this loop allows you to visualise multiple plots. In this case, it will run through the first 4 spectra and peaks, and create two plots for each, one entire spectra, and one zoomed in to a region of interest
for (i in 1:4){
  plot(spectra[[i]])
  points(peaks[[i]], col="red")
  plot(spectra[[i]], xlim=c(1500, 1600)) #xlim sets the limits of the x-axis. Edit to zoom in on regions of interest
  points(peaks[[i]], col="red")
}
```
## Remove low-quality spectra from technical replicates

```{r}
x <- c() #empty vector for number of peaks
n <- length(peaks) #number of spectra analysed

for (i in 1:n){
  x[i] <- length(peaks[[i]]@mass) #creates a list with count of peaks for each spectra
}

hist(x, breaks = 30) #View shape and distribution of the data including the number of peaks and modes

remove1 <- which(x<45) # index location of peaks less than 45 (depending on histogram)
remove2 <- which(x>105) # index location of peaks less than 105 (depending on histogram)

remove <- c(remove1, remove2) #creates an object with index locations to remove from spectra sample set
```
# Creates subset of high-quality spectra -- merging the high quality technical replicates
```{r}
spectraHQ <- spectra[-remove] #this removes individual spectra in the object called "remove"
```
# Averaging replicates

## In our files: "ID#_spot.mzML" -- run this following code snippet if you have the same format, otherwise adjust in code

```{r}
n<-length(spectraHQ)

for (i in 1:n){
  file <- basename(spectraHQ[[i]]@metaData[["file"]])
  ID <- gsub("_.*", "", file) #this removes everything after the "_" for averaging
  spectraHQ[[i]]@metaData[["ID"]] <- ID #creates a new line of metadata with sample IDs
}
```
## Create a list of the ID names
```{r}
grouplist <- c() #create empty vector object

for (i in 1:(length(spectraHQ))) {
  grouplist[i] <- as.character((spectraHQ[[i]]@metaData[["ID"]]))
}

grouplist <- as.factor(grouplist) #convert to a factor
```
```{r}
spectraHQ<-averageMassSpectra(spectraHQ, labels=grouplist, method="mean") #averaging high technical replicates
```
## Peak pick and deisotope the high quality averaged spectra

```{r}
peaksHQ <- detectPeaks(spectraHQ, method="SuperSmoother", halfWindowSize=20, SNR=3.0) #change SNR
peaksHQ <- monoisotopicPeaks(peaksHQ, minCor=0.95, tolerance=1e-4, distance=1.00235, size=2L:10L) #adjust parameter values based on your spectra -- see README.md in SpecieScan R repo

## Visualise peaks
for (i in 1:10){ 
  plot(peaksHQ[[i]])
  plot(peaksHQ[[i]], xlim=c(1500, 1600)) 
}
### If error out of bounds -- go back to the beginning of the code and change parameter values or make the range of remove1 and remove2 larger
# Export spectra as csv (create folder for it first)

exportCsv(peaksHQ, path= "~/Documents/Site_name_preprocessed_csvs") #you move on to the SpecieScan Python code with these files only 
