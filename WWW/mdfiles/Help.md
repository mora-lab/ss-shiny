### **1. What is SS-shiny?**
#### SS-shiny is a shiny app with the goal of benchmarking different single-sample Gene Set Analysis methods in an easy way.
#### Among its main advantages, SS-shiny allows the user to: (i) re-define the information that is using as a gold-standard (test different gold standards), 
#### (ii) generate results of the ensemble of methods (through combining individual p-values), 
#### (iii) Choose the benchmarking metrics of interest (sensitivity, specificity, precision, ROC).

### **2. What is the correct format to upload RNA data?**
#### The user can work with a sample data or the gold standard built by Tarca et al., 2013, which are already in the correct format. 
#### If the user chooses to upload its own gold standard dataset, the data must be uploaded as a single .RDS file (created in R). 
#### The .RDS file should contain a list of expression data matrices (each element of the list being a different RNA dataset). 
#### Each expression data matrix should include rows representing gene symbols (except for the first row) and columns representing samples. 
#### The first row should contain the information on each sample's disease status, that is, either 0 (disease/tumor) or 1 (normal/control). 
#### You can check our sample dataset for more information.

### **3. What is the correct format to upload pathway data?**
#### Same as with the previous question, the user can work with our sample file or Tarca's list of target pathways, which are already in the correct format.
#### If the user wants to upload its own pathway information, the format should also be an .RDS file, this time containing a list of pathways, each pathway 
#### made of a vector of gene symbols. You can check our sample dataset for more information.

### **4. Why only 5 GSA mehods? What if I want to compare the results of a different method?**
#### SS-shiny only includes methods for the SS-GSA category which are also reasonably fast. Testing additional methods is not possible through the shiny interface, 
#### but they are possible to be added using the underlying R infrastructure. To see an example of how can that be done, you can read the following jupyter notebooks:
#### * [GSVA workflow](https://github.com/mora-lab/benchmarks/blob/master/single-sample/workflows/GSVA_GSE10245.ipynb) 
#### * [Pathifier workflow](https://github.com/mora-lab/benchmarks/blob/master/single-sample/workflows/pathifier_GSE10245.ipynb)
#### * [GRAPE workflow](https://github.com/mora-lab/benchmarks/blob/master/single-sample/workflows/GRAPE_GSE10245.ipynb)
	
### **5. Why should I choose a method to combine p-values?**
#### As it was said before, SS-shiny will also generate results from an ensemble of all methods included in the analysis.
#### This is done by combining the resulting p-values using one of the specified statistical methods.

### **6. What are sensitivity, specificity, precision, and ROC plot?**
#### Sensitivity measures the proportion of true positives that are correctly identified as true positives, TPR = TP/(TP+FN). 
#### Specificity measures the proportion of true negatives that are correctly identified as true negatives, TNR = TN/(TN+FP).
#### Precision is used to evaluate the quality of a classification system, answer the question:"How close the values that have been taken are close to each other?", precision = TP/(TP+FP).
#### ROC stands for Receiver Operating Characteristic (from Signal Detection Theory), ROC Curves can be used to evaluate the tradeoff between true- and false-positive rates of classification algorithms.
#### here it needs only TPR(Sencitivity) and FPR(1-Specificity).

### **7. How to download the results?**
#### In the "Data analysis" section, you must start by introducing all required information, as explained before. Then you will have access to three different tabs: 
#### "Preview of datasets" will show you a preview of the input; 
#### "Preview of results" will show you the resulting comparison plots; 
#### finally, "Downloads" will allow you to download such results.

### **8. How was SS-shiny built?**
#### We have developed multiple R functions and datasets for benchmarking Gene Set Analysis methods. That information can be found here: https://github.com/mora-lab/benchmarks/tree/master/single-sample
#### SS-shiny is a web-style shiny interface to those functions which is limited to five SS-GSA methods. The open code can be found here: https://github.com/mora-lab/benchmarks/tree/master/single-sample/SS-shiny
#### For alternative ways to use our benchmarking functions, find some benchmarking workflows (as jupyter notebooks) here: https://github.com/mora-lab/benchmarks/tree/master/single-sample/workflows

### **9. If facing problems, who should I contact?**
#### For questions about specific applications or the general scope of the project, contact Antonio Mora (antoniocmora@gzhmu.edu.cn)
#### For questions regarding technical problems or specific bugs, please contact Chengshu Xie (chengshu@std.gzhmu.edu.cn)

 
	