# Thylakoid_VD
Python Code to analyze thylakoid grana stack vertical dimensions from x,y data of TEM images

####====Software requirement====####

Python 3.7, Text editor Spyder 4.0.0 or above recommended
The code has been tested last time in Python 3.7. It should also work in any Python 3 environments.
Besides, ImageJ or Fiji software was used to generate "x,y data" from TEM images.

####====Input and Output====####

ThylakoidVD_csv.py was written to analyze the vertical dimension of thylakoid grana stacks. With 
input of one or more X,Y csv files (tsv file should work too), for each input csv file: output one txt file of 
calculated Repeat distances, Stroma gap width... (* LSRD.txt) and the * plots.pdf file of data and fit to inspect the quality
of fittings. One should NOT use the generated LSRD.txt file as it is to do statistical analysis, since usually 
some errors would occur either on the curve fitting level.

####====recommended way of using the code====####

After compiling the x,y data files (in cvs format), place those csv files together with ThylakoidVD_csv.py in the
same folder and run the ThylakoidVD_csv.py within the folder. It should read all the csv files and try to prossess
all the files. 
One can download the Test_ThylakoidVD_csv.zip file to test run the code with some real data generated from our HPF study.
Test_ThylakoidVD_csv.zip folder has a copy of the Python code and experimental data for test. Running the code will generate
the pdf and txt files given in the same folder. One can delete or move them (txt and pdf files) out to see if it works.

####====The basic design====####

The calculation principle is described in the supplemental method 1 section of the associated manuscript:
----------Not yet published, to be updated---------

####====How to get the correct input csv file?====####

Just in case one has not dealt with similar problem or never studied thylakoid structure, please see the following illustrated 
screenshot for how to generate an input "x,y" dataset using ImageJ or Fiji:

How_to_prepare_for_ThylakoidVD_csv.png

