
Here we provide code used in the following papier: 
"Penetrance estimation of Alzheimer disease in SORL1 loss-of-function variant carriers using a family-based strategy stratified by APOE genotypes"
By Catherine Schramm, Camille Charbonnier, Aline Zaréa, Morgane Lacour, David Wallon, CNRMAJ collaborators, Anne Boland, Jean-François Deleuze, Robert Olaso, Flora Alarcon, Dominique Campion, Grégory Nuel, Gaël Nicolas
doi: https://doi.org/10.1101/2021.06.30.450554 

 example.txt is a simulated dataset  
 
 example.R is the R file able to estimate the age-related penetrance of Alzheimer disease at the digenic level (SORL1 and APOE)
 
During the estimation process, the R function needs to call a C++ implementation of bped performing the E part of the EM algorithm (posterior genotype probability estimation). To be able to use the entire function, the user should first download the bnlib2.h and bped_smallvar_3alleles_2alleles.cpp files and compile them following instruction provided by these files. Then the user should modify the line 145 of the "run_model.R" file to provide the path to the bped3alleles2alleles C++ function obtained after compiling the bped_smallvar_3alleles_2alleles.cpp file.

For an easier use of our method, an R version of bped will be available soon.

### Additional information

For questions about the method, please contact us: cath.schramm@gmail.com

For bped use, please note the following license:

copyright G. Nuel (2014-2022)  
This program is distributed in the hope that it will be useful, but without any warranty. It is free to use for academic research, but not for commercial or medical purposes. For any question about this code, please contact the author (nuel@math.cnrs.fr).
