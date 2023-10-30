# compile

     make

When the 'make' command is finished,'FeatureGen' and 'Cluster' files are generated.    

# config file 
  The parameters for the procedure are recorded in config file. 
  
The content of config file:
    
    Tolerance = 0.000001
    CutOff = 0.20
    MaxDepth =  1024
    MinSupport = 3
    minRatio = 0.1
    VarianceThreshold = 1.0

Tolerance:

   ![first equation](http://latex.codecogs.com/gif.latex?%5Csigma%20_%7Bi%7D) is an eigen value of matrix.  
   
   If Tolerance is specified, threshold = 1-Tolerance; The procedure will find the maximum r which satisfy  
   ![first equation](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Br%7D%5E%7B2%7D%7D%7B%5Csigma%20_%7B1%7D%5E%7B2%7D&plus;%5Csigma%20_%7B2%7D%5E%7B2%7D&plus;%5Csigma%20_%7B3%7D%5E%7B2%7D%20...%20&plus;%5Csigma%20_%7Bn%7D%5E%7B2%7D%7D%5Cleq%20threshold). Then the first r principal components are selected. 
 
CutOff: 
   
   There is an F-distribution for all attribute values in a component. The P-value of each selected attribute calculated by F-distribution must be smaller than CutOff. So the number of selected attributes in each component depends on CutOff. Smaller CutOff can shrink the number of final results. If the number of results exponentially expands too large, the CutOff should be tuned smaller.
   
MaxDepth:

   The maximum number of principal components of Attribute Clustering. 

MinSupport:

   The minimum size of Attribute Clustering. The bigger MinSupport will introduce a smaller number of final results.
   
minRatio: 
   
   The ratio is the percentage of the selected components for a cluster to all components for an attribute. This is the minimum ratio for all attributes of a cluster.

VarianceThreshold:

    If the variance for a component divided by the singular value of this component is smaller than the Threshold, then this component is removed. This threshold is used to remove the little change components. 

# Comands

    ./Cluster TestMatrix   OutFile   [ fileName ] 
    
The config file is necessary that exist in the same directory as the data file. if not, the path of the config file must be specified.

    ./Cluster TestMatrix  OutFile  [ fileName ]  -config   configPATH/config 
    
'fileName' is optional, the column names of the TestMatrix file are stored in this file, the format of the 'fileName' is described as follows: 

    column_name1
    column_name2
    ...
    column_namen
    
 If this file is ignored, the names of the TestMatrix file are as named as '0','2','3'..'N', N is the length of 'fileName'. 
         
'TestMatrix' is the input file, the format of the input file is described as follows: 

      0   0   4.00
      0   1   3.00
      0   2   4.00
      0   4   7.00
      0   5   111.00
        ... 
        
 The first column indicates the column position and the second column indicates the row position. The third column indicates the value of this position. 


 'OutFile' is output file name. 
 
        grep CorrelationAverage OutFile > log 
  
 After the above command is executed, the content of the log file is described as follows: 
 
        CorrelationAverage: 3 0.999124 0.148882
        CorrelationAverage: 3 0.989092 0.142556
        ......
  The first column numeric indicates the number of elements in each clustering. The second column indicates the Average correlation of this clustering and the third column indicates the Average Component Ratio of this clustering.  
  The  distribution of Average correlation and Average Component Ratio can be obtained from the log file, The thresholds of Average correlation and Average Component Ratio for filtration are evaluated from this distribution. 
 
      
     FeatureGen ToFeature   InputFile OutputFile AverageCorrelation AverageComponentRatio
 
 'ToFeature' is the keyword. InputFile is the filename of the Cluster command-generated file.  OutputFile is the output file name for the features record. AverageCorrelation is a numeric to indicate the threshold of Average correlation and AverageComponentRatio is a  numeric to indicate the threshold of Average Component Ratio.  
 
    FeatureGen Decipher FeatureFile   InputFile OutputFile  [-N TableNameFile ]
   
  'Decipher' is the keyword, 'FeatureFile' is the filename of 'FeatureGen ToFeature' command generated file. 'InputFile' and 'OutputFile' are the matrix files with the same format as the input file ('TestMatrix') of 'Cluster' command.  'InputFile' is the input matrix and  '-N TableNameFile' indicates the 
column name of the input matrix. The format of TableNameFile is the same as 'fileName' file of 'Cluster' command. 
  
  

