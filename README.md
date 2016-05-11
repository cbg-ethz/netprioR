# netprioR

*netprioR* is a probabilistic graphical model for gene prioritisation. The model integrates network data, such as protein--protein interaction networks or co-expression networks, gene phenotypes, e.g. from perturbation experiments, as well as prior knowledge about a priori known true positives and true negatives for a prioritisation task. The goal of the model is to provide a robust prioritisation (i.e. ranked list) of genes accounting for all dependencies in the input data. Parameter inference and imputation of hidden data is performed in an Expectation Maximisation (EM) framework.

Check out the vignette for a showcase of the main functionality.
