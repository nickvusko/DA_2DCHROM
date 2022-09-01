# DA_2DCHROM
Software for alignment of two dimensional chromatograms

Readme for DA_2DCHROM

	####CONTACT INFORMATION####

	Author: 	Nikola Ladislavová	(ladislan@vscht.cz; ORCID 0000-0001-8733-4780)
				UCT Prague, Department of Analytical Chemistry, UCT Prague

#########################USER GUIDE#########################
VERSION: 1.0

REQUIREMENTS
Python == 3.9.7

Matplotlib == 3.4.3
Numpy == 1.20.3
Pandas == 1.3.4
Scipy == 1.7.1



	####OVERALL SCRIPT DESCRIPTION####

DA_2DCHROM.py is designed for the needs of two-dimensional chromatographic data manipulation (data alignment). The input data structure corresponds to the data structure of LECO ChromaToF(R) .csv export format.
DA_2DCHROM.py offers 5 different data alignment approaches - BiPACE2D, DISCO, MSort, PAM, RI, and TNT-DA methods along with adjustable method for comparing spectral similarities. The script should be used for automatic data alignment of multiple chromatograms.


	####INPUT DATA FORMAT####
The script takes a .txt file with a tabulator as a separator. The files must contain columns*: "1st Retention Time(s)", "2nd Retention Time(s)", and "Spectra" or "Spectrum". The script does not drop or edit any additional columns; the listed columns are essential for data alignment process. "1st Retention Time(s)" should contain integers (or floats), "2nd Retention Time(s)" should contain floats, and "Spectra"/"Spectrum" should be in following string format: mass:intensity (blank_space) mass:intesity.

Sample of the desired file format:
1st Dimension Time (s)	2nd Dimension Time (s)	Spectra
668	2.5	73:999 281:683 45:358 75:343
710	3.175	269:998 29:872 267:739 31:631
764	3.73	43:999 57:748 71:627 41:620
794	2.87	73:999 341:144 45:133 75:125

*If you are importing files with different column names, you can always rename the columns by the following command:
try:
	file = file.rename(columns=“column_name_to_rename”: “desired_column_name”)
except KeyError:
        print("Variable names are already correct .")
        pass

Please note that you need to put the names between the quotation marks as shown in the example. Insert the lines shown to line nmr. 280 (function merge_peak_entries), just after the block of renaming the "Spectra" column.

	####EXECUTION PATTERN OF DATA ALIGNMENT####

The general path of the script is as follows:

Find anchor peaks -> Check elution order of anchor peaks in both dimensions -> Export anchor points -> Retention time shift correction -> Export peak map -> Export time corrected chromatogram

The “Find anchor peaks" part depends on the data alignment method and spectral comparison method chosen by the user.


	####ARGUMENTS OF THE DA_2DCHROMSCRIPTS####

The main "run" line 964 of the script:
x=data_align(Ref_chrom :str, Preprocess=True, Source_folder = "data", Result_folder= "metadata/time_correction", method="DISCO", similarity="Pearson", transform =(1.,0.)).run_process()

##Ref_chrom## -> str
The only mandatory input from the user is the Ref_chrom argument. It needs to be in a string format ("name_of_the_file_without_a_file_extension").
For example, if the line looks like this:
data_align("500_system_1_M13_4")
The script will take that file as a referential chromatogram, performs pre-processing step, it takes all the files in a "data" folder, saves the result into the "metadata/time_correction", the alignment method will be "DISCO", the similarity of spectra will be calculated by Pearson correlation coefficient, and the spectral data will not be transformed in any way.
Please note that the Ref_chrom file needs to be in the Source_folder.

##Preprocess## -> bool
The preprocess argument considers True/False values. It is set to True by default. The preprocessing step takes the originial data, checks the number formating (replaces "," for "." in numerical columns), renames columns into the proper name format, and merges peaks of the same chemical substances. The preprocessed data are stored in the "metadata/merged_peaks" folder, and the data alignment itself is performed on the preprocessed data. If you set pre-process to False, the script will perform the data alignment straight on the files stored in Source_folder. The Pre-process step also deals with the ChromaToF feature which adds a "saturated (number )" comment into the "Area" and "Height" column.

##Source_folder## -> str
The Source_folder argument should be a full existing folder path (recommended). If you are getting the "File not found" error, please check the folder path. The default value is set to working_directory**/data folder.

##Result_folder## -> str
The Result_folder argument should be a full folder path (recommended). If the folder does not exist, the script will create the folder.

** You can check your working directory by following commands:
import os
os.getcwd()

##method## -> str
For version 1.0, there are 6 methods available: BiPACE2D[1], DISCO[2], MSort[3], PAM[4], RI, and TNT-DA
See the references for the explanation and selection rules of the parameters of each algorithm.
	#BiPACE2D#
The algorithm takes 5 arguments: RT1 tolerance parameter, RT2 tolerance parameter, T1 critical value, T2 critical value and similarity threshold value.
Setting of T1, T2 influences the decision whether the algorithm should proceed with peak comparison or whether the two compared peaks are too far away.

	#DISCO#
The algorithm takes 1 argument: similarity threshold value.
To speed up the processing speed, the algorithm considers only 20% of the total rows with the least distance for each ref_table row.
This parameter can be changed by the bound = int((len(align_chrom)/100)*20) argument (line 440)

	#MSort#
The parameters max_shift_1 and max_shift_2 are the max allowed deviation of peaks in the 1st and 2nd dimension (please note that the scale of 1st and 2nd differs). max_shift_1 should be a multiple of the modulation period(s)
The similarity parameter is the threshold for the mass spectra similarity of compared peaks.

	#PAM#
The parameter w is crucial for the calculation of mixture similarity: w / (1 + distance) * (1-w) * similarity.

	#RI#
The RI method is built on manual identification of anchor points. The algorithm is specifically designed for measuring and evaluating scent samples at UCT Prague and it is not recommended to use this method for data alignment in general. 
RI has a rigid input format as shown below:
Name	R.T. (s)
Ref1_1	704 , 3.490
Ref1_2	974 , 3.630
Ref1_3	1220 , 3.750
Ref1_4	1442 , 3.810
Ref1_5	1634 , 3.890
Ref2_1	1676 , 1.190
.......

First, anchor points need to be exported from ChromaToF(R) software: 2 columns - "Name" and "R.T. (s)". There must be 17 anchor points with specific names in total - 13 for 1st dimension (Ref1_1....Ref1_13), and 4 for 2nd dimension 
(Ref2_1, Ref2_2, Ref2_3, Ref2_5 - Ref2_4 refers to Ref1_8 in the specific scent sample setup). All files containing reference data should be stored in "metadata/ref_data/" folder to execute the scrip correctly. The script simply calculates Retention Indexes of all peaks in the chromatogram. Please note that it is only an approximation; thus, technically, it is not a proper data alignment method.

	#TNT-DA#
TNT-DA is derived from DISCO discrimination rules and the PAM selection rule. First, retention times are recalculated to z-scores. For every peak in a reference table, Canberran (or Euclidian) distances for every peak in an aligned table are calculated. Then 20 % of the nearest peaks are considered as potential anchor peaks. Spectral similarity (cosine similarity with 0.53, 1.3 transformation) is calculated for every potential anchor peak. If the spectral similarity overreaches the 90 % (the value is set by the user) similarity threshold, PAM selection rule is applied, and the peak from aligned table with the highest score is marked as an anchor peak.

##similarity## -> str
There are two available methods for spectral similarity comparison - Pearson (Pearsons correlation coefficient) and DOT (cosine similarity).

##transform## -> float
transform =(a,b) takes the transformation coefficients for the mass spectra. The intensity for each mass is recalculated by the following formula:
((intensity)**a)*(m/z)**b

	####OUTPUT DATA FORMAT####
Every aligned file produces 3 files – aligned chromatogram (file name format: method_parameters_originalfilename.txt), list of anchor peaks (file name format: file name format: method_parameters_originalfilename_anchors.txt), and additional shift map (file name format: file name format: method_parameters_originalfilename.png)
The aligned file format is identical to the input file format. Anchor peaks file format refers to reference peak ID (very first column) along with its time coordinates (“1st RT_ref” and “2nd RT_ref” columns) and aligned peak ID (“match”) along with its time coordinates (“1st Dimension Time (s)” and “2nd Dimension Time (s)” columns)

	####RUNNING THE SCRIPT####
You can run the code simply through the cmd: python “path_to_the_script/DA_2DCHROM.py”.
If you want to edit arguments, you need to open the code in any .txt editor (or other software) and simply rewrite the argument.
	####LICENSE####
This script is licensed under Creative Commons Attribution 4.0 International.

####ACKNOWLEDGMENT####
The software development was funded by TACR (programme IMPAKT 1; project VJ01010123)
	####REFERENCES####
[1]	Hoffmann, N., et al., BiPACE 2D—graph-based multiple alignment for comprehensive 2D gas chromatography-mass spectrometry. Bioinformatics, 2013. 30(7): p. 988-995.
      DOI: https://doi.org/10.1093/bioinformatics/btt738

[2]	Wang, B., et al., DISCO: Distance and Spectrum Correlation Optimization Alignment for Two-Dimensional Gas Chromatography Time-of-Flight Mass Spectrometry-Based Metabolomics. Analytical Chemistry, 2010. 82(12): p. 5069-5081.
      DOI: 10.1021/ac100064b

[3]	Oh, Cheolhwan, et al. "Comprehensive two-dimensional gas chromatography/time-of-flight mass spectrometry peak sorting algorithm." Journal of chromatography A 1179.2 (2008): 205-215.
      DOI: https://doi.org/10.1016/j.chroma.2007.11.101

[4]	Kim, Seongho, et al. "An optimal peak alignment for comprehensive two-dimensional gas chromatography mass spectrometry using mixture similarity measure." Bioinformatics 27.12 (2011): 1660-1666.
      DOI: https://doi.org/10.1093/bioinformatics/btr188
[5] doplnit
#########################PROGRAMMER GUIDE#########################
####Class data_align()####
	Initial attributes are set to the chosen methods accordingly. The script will prompt the user to type in the necessary information.

##__init__(self,Ref_chrom :str, Preprocess=True, Source_folder = "data", Result_folder= "metadata/time_correction", method="DISCO", similarity="Pearson", transform =(1.,0.), multipool=True)##
	Initialize the attributes of the class objects. 

	##run_process (self)##
run_process function is the core function of the script. It initializes every step of the data alignment process. First, preprocessing bool is evaluated and executed. Then, both referential and aligned tables are loaded, compared, and the elution order in both dimensions is checked. Finally, anchor points, shift maps, and aligned chromatograms are exported to files described in output section.

	##data_align_chrom (self, lst_element:str)##
The function loads compared chromatogram and finds the anchor points. Then, it checks their elution order in both dimensions, exports the list of the anchor peaks, shifts the compared chromatograms and exports aligned chromatograms with the alignment maps.

	## compare_chroms (self,align_chrom, method) ->pd.DataFrame## 
compare_chroms() function is signpost. The function identifies the align method and executes the corresponding script. The function returns the list of found anchor peaks.

	##merge_peak_entries (self, filename:str)-> pd.DataFrame##
This is a pre-process function. First, it checks if whenever the aligned file is in proper format. The lines transform referential peak spectrum to a dictionary along with the spectral intensity transformation. Subsequently, time search window is defined (with respect to the retention drift), and peaks found in the search window are compared based on their mass spectra. If the similarity is higher than set threshold, peaks are marked as the same chemical compound and merged. The function exports merged dataframe.

	## compare_spectra (self, samplemasses_main:dict, spectrum2:list)->float##
compare_spectra () transform the spectra of an investigated peak in the same way as lines. Afterwards, both dicts are put into a single dataframe and compared via Pearson or DOT method. The function returns a float result in %.

	## mask_groupby (self,df:pd.DataFrame)->pd.DataFrame##
This is an inner function of merge_peak_entries (). After the outer function “tags” peak of the same chemical compound via indexes, mask_groupby () function simply merges rows with the same “index” value. The area and height columns are summed, and retention time coordinates are recalculated. The new merged peak spectrum is set as a spectrum from the original row with the highest area value. The function returns a new merged dataframe.

	##calculate_zscore_times(self,dataframe:pd.DataFrame)-> pd.DataFrame##
This function is triggered by alignment methods requiring z-score transformation. The function adds two new columns to the dataframe: "Z score 1st RT" and "Z score 2nd RT". The function returns a new dataframe

	##check_rank (self, anchor_points:pd.DataFrame, ref_table:pd.DataFrame)-> pd.DataFrame##
This function checks whether an elution order of the peaks in the first dimension is same in both referential chromatogram and aligned chromatogram. If there is a difference in the elution order, anchor peak is removed as it is the result of the wrong spectral identification. The function returns a new dataframe containing the list of anchor peaks.

	##shift_1stRT(self,anchor_points:pd.DataFrame,ref_chrom:pd.DataFrame, align_chrom:pd.DataFrame)->pd.DataFrame##
The function takes the anchor peaks dataframe, aligned chromatogram dataframe, and alignes all data points in the aligned chromatogram according to the shift of the anchor peaks through the linear approximation. shift_2ndRT is analogous function to shift_1stRT for the second dimension. The function exports aligned chromatogram dataframe.
