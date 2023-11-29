## You can install the package via:

	library(devtools)
	install_github("NovakLab/FL_Bstate")

## Usage
The package requires normalized gene expression values. The methodology will automatically detect if a log transformation has been applied and adjusts accordingly. The gene expression values can be from light based or sequencing based technologies. The signature matrix is derived from RNAseq TPM values, so batch correction is advised when using other technologies. The row names need to be either Hugo symbols or ensembl gene IDs. If you do not have those values, you will need to convert your input to have them.  

The data is contained in:
  
	data("test_exp") #test FL expression
	data("signature_mat") #signature matrix for prediction  
Use the `signature_mat` to make sure your gene expression matrix is in the correct orientation, but it should be genes as rows and samples as columns.  

 	

Basic functionality:  

	test_pred_df <- Coef_predictor(test_exp, nu = 0.1, batch_correct = FALSE)  

The result of this function is a dataframe of samples as rows and columns `c("INFM", "PDZ", "CMI")`. The values are scaled state values for the three states and each sample.  
	
It isn't recommended to mess with the nu value unless you are familiar with SVR and what changing the nu value means. The default value is set to the optimum for this study. As was mentioned, the batch correction is to be set to TRUE for non-RNA-Seq data.

### Class Prediction
You can also use an additional function to predict the class label. Once you've predicted the state values, you can use this function to predict the dominant signature.

	predict_class(temp_df)

The result of this is a named list of class labels for each sample.

[Link to HTML file](https://github.com/NovakLab/FL_Bstate/blob/main/docs/index.html)
