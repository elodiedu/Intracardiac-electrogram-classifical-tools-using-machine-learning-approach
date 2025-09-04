The matlab files shows how to extract the feature data, MSF, MSE .etc from the signals exported from each stage; The imported data used are from the Matfile folder. Custom algorithms need to be written to get the format like the data in the Matfile folder from the Raw Carto 3 data. 

Run the file convert mesh to visualize the 3D mesh. The iEGMsAnalysis is the main file you should run, as inside it calls the other function such as NMSE NNMA and MI, keep them in the same folder. The iEGMsAnalysis should export extracted features, as in the format of the data folder. There are two types of format, the iEGMsinformation.. which compile the three features in the same matrix, and the RA1MSF, which contains only one feature per file. Select which data you want to use based on the need, the difference is only the format.


Now open the .py file to do classifier. The Classifier.py is the main file, which shows multiple machine learning classifier and there visualization result. The K_fold test is a supplementary file that provide additional information about the K-fold test result. 
