# Integrating LASSO feature selection and soft voting classifier to identify origins of replication sites
DNA replication plays an indispensable role in the transmission of genetic information. It is considered to be the basis of biological inheritance and the most fundamental processes in all biological life. Considering DNA replication initiates with a special location namely the origin of replication, the better and accurate prediction for the origins of replication sites (ORIs) is essential to gain insight into the relationship with gene expression. In this study, we develop an efficient predictor called iORI-LAVT for ORIs identification. This work focuses on extracting feature information from three aspects including mono-nucleotide encoding, k-mer and ring-function-hydrogen-chemical properties. Subsequently, least absolute shrinkage and selection operator (LASSO) as a feature selection is applied to select the optimal features. Comparing the different combined soft voting classifiers results, the soft voting classifier based on GaussianNB and Logistic Regression is employed as final classifier.


#DATASET
There are three datasets, where S1  and S2 are training datasets and S3 is an independent dataset. The dataset S1  consists of 405 ORI sequences and 406 non-ORI sequences, which belongs to Saccharomyces cerevisiae (S.cerevisiae) genome. S2  and  S3 are derived from Arabidopsis thaliana (A. thaliana) genome, where S2 is made of 1015 ORI sequences and 1015 non-ORI sequences,  S3 is made of 500 ORI sequence and 500 non-ORI sequences. Moreover, the length of all given DNA sequences is 300bp.
