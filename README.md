# SASNE - Shape-aware Stochastic Neighbor Embedding


### Introduction 

SASNE (Shape-aware stochastic neighbor embedding) is a method that can be used to create low dimensional visualisation that reveal key structures of high dimensional data.

For more information about SASNE and possible applications please refer to our [preprint](https://doi.org/10.21203/rs.3.rs-1831618/v1). 

The method is currently implemented in MATLAB.

### Installation

To install SASNE input the following command to a terminal:

    git clone  https://github.com/tobiaswangberg/SASNE.git
    cd SASNE
    
Alternatively download the project as a ZIP file directly from the `github` page.

### System requirements

The program requires a [MATLAB](https://www.mathworks.com/products/matlab.html) installation (version >= 2019a).

### Tutorial

The SASNE function takes as input a matrix, where rows are observations and columns are features. The program returns the two dimensional embedding. 

For an example of how to apply SASNE run the following code in your MATLAB instance:

    data = readmatrix('data/imbalanced_test.txt'); 
    [sasne_out,Z] = SASNE(data);

To visualise the low dimensional embedding run the following command in MATLAB:

    scatter(sasne_out(:,1),sasne_out(:,2),5,'filled')
    
Since not all data can be accurately embedded into 2 dimensions it is important to include a validation measure to quantifies to what extent the original structures in the high dimensional data is preserved in the low dimensional embedding. To this end include a Rank Residual Plot (see [our paper](https://doi.org/10.21203/rs.3.rs-1831618/v1) for details) using the following code

    D_original = pdist(Z,'euclidean');
    D_sasne = pdist(sasne_out,'euclidean');
    RRP(D_original,D_sasne);

The same running example is found in the live script `running_example.mlx` found in this repository.

### Help

If you have any general questions or need help with running the program please feel free to send a message to <tobias@math.su.se>.
