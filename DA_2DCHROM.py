# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:06:13 2022

@author: ladislan
"""

import pandas as pd
from scipy import stats
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import math
from numpy import dot
from numpy.linalg import norm
import concurrent.futures
import sys

from warnings import filterwarnings
filterwarnings("ignore")

class data_align():
    """
    Data_align class object for processing of chromatograms via various methods.
    Containts all basic processes and preprocesses such as merging same peaks, mass spectra compare algorhitms and core structure for data alignemnt.
    
    
    The initial parameters for the data_align class:
        Ref_chrom       : str
                        The name of a csv file which serves as a reference talbe
                        "example_of_ref"
        
        Preprocess      : bool
                        If False, the script will skip preprocessing steps (formating the files, merging multiple peaks of the same chemical origin through the processed chromatogram)
                        Default value = True
        
        Source_folder   : str
                        Path to the input chromatograms
                        Default value = "data/"
        
        Result_folder   : str
                        Path to the output chromatograms
                        Default value = "metadata/time_correction"
                        
        method          : str
                        Data alignment method
                        Default value = "DISCO"
                        Available methods: BiPACE2D, DISCO, MSort, PAM, RI, TNTDA
                        DISCO: Similarity treshold
                        BiPACE2D: D1, D2, T1, T2
                        MSort: 1st dimension max shift, 2nd dimension max shift, similarity treshold
                        PAM: The weight for the mixture similarity measure, method for calculating distances between compared peaks.
                        TNTDA: The weight for the mixture similarity measure, similarity treshold, method for calculating distances between compared peaks                
        similarity     : str
                       Method for calculating spectral similarity
                       Default value = "Pearson"
                       Available methods: Pearson, DOT
        
        transform      : tuple(float(a),float(b))
                       The transformation equation of mass spectra: (intensity)**a*(m/z)**b
                       Recommended values are (0.53,1.3) *Kim, Seongho, et al. "Compound identification using partial and semipartial correlations for gas chromatography–mass spectrometry data." Analytical chemistry 84.15 (2012): 6477-6487.
                                                      DOI:https://doi.org/10.1021/ac301350n
        
        multipool     : bool
                      Use multiprocessing approach
    """
    
    def __init__(self,Ref_chrom :str, Preprocess=True, Source_folder = "data", Result_folder= "metadata/time_correction", method="DISCO", similarity="Pearson", transform =(1.,0.), multipool=True):
        if not os.path.isdir(Result_folder):
            os.makedirs(Result_folder)
            print(f"The result folder {Result_folder} was created.")
        self.Lst_of_files = [file for file in os.listdir(f"{Source_folder}") if file.endswith(".txt")]
        self.Result_folder = Result_folder
        self.Ref_chrom_name = Ref_chrom
        self.Preprocess = Preprocess
        self.Source_folder = Source_folder
        self.Method = method
        self.Similarity = similarity
        self.Transformation=transform
        self.multipool = multipool
        if self.Preprocess == True:
            self.minimal_merge = 60 #int(input("What is the minimum similarity (in %) for peak merge matching?"))
            self.Folder_merged = f"{self.Result_folder}{os.sep}merged_peaks"
            if not os.path.isdir(self.Folder_merged):
                os.makedirs(self.Folder_merged)
        if method == "DISCO":
            self.similarity = 60 #float(input("Similarity treshold parameter: "))
        elif method == "BiPACE2D":
            self.D1 = float(input("D1 parameter: "))
            self.D2 = float(input("D2 parameter: "))
            self.T1 = float(input("T1 parameter: "))
            self.T2 = float(input("T2 parameter: "))
            self.similarity = float(input("Similarity treshold parameter: "))
            print(f"D1: {self.D1}, D2: {self.D2}, SI: {self.similarity}")
        elif method == "MSort":
            self.max_shift_1 = float(input("1st dimension retention time max shift parameter: "))
            self.max_shift_2 = float(input("2nd dimension retention time max shift parameter: "))
            self.similarity = float(input("Similarity treshold parameter: "))
        elif method == "PAM":
            self.w = float(input("w parameter: "))
            self.distance = "Canberra"#str(input("Distance measure, choose one: [Euklidian, Canberra]"))
        elif self.Method == "RI":
            pass
        elif self.Method == "TNTDA":
            self.w = float(input("w parameter: "))
            self.similarity = float(input("similarity threshold parameter: "))
            self.distance = str(input("Distance measure, choose one: [Euklidian, Canberra]"))
        else:
            print("No valid method was given as an input. Program will be terminated")
            quit()
            
    def run_process(self):
        """
        The sequence of functions for finding anchor points, relaclculation of time shifts and export of the results with both numerical and graphical representations. 

        Returns
        -------
        align_chrom : pd.DataFrame
            The final aligned chromatogram

        """
        if self.Preprocess == True:
            start = datetime.now()
            if self.multipool == True:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = [executor.submit(self.merge_peak_entries, chromatogram) for chromatogram in self.Lst_of_files]
                    for f in concurrent.futures.as_completed(results):
                        print(f.result())
            else:
                for chromatogram in self.Lst_of_files:
                    self.merge_peak_entries(chromatogram)
            print("Execution time:", datetime.now()-start)
        if self.Preprocess == True:
            self.Lst_of_files = [file for file in os.listdir(f"{self.Folder_merged}") if file.endswith(".txt")]
            self.Source_folder = self.Folder_merged
        self.Ref_chrom = pd.read_csv((f"{self.Source_folder}{os.sep}{self.Ref_chrom_name}.txt"), sep = "\t", header = 0, encoding= "latin-1")
        if self.multipool == True:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = [executor.submit(self.data_align_chrom, lst_element) for lst_element in self.Lst_of_files if lst_element != self.Ref_chrom_name and self.Method != "RI"]
                for f in concurrent.futures.as_completed(results):
                    print(f.result())
        else:
            for lst_element in self.Lst_of_files:
                self.data_align_chrom(lst_element)
            
    def data_align_chrom(self, lst_element:str):
        start_time = datetime.now()
        self.align_chrom_name=lst_element
        align_chrom = pd.read_csv((f"{self.Source_folder}{os.sep}{lst_element}"), sep = "\t", header = 0, encoding= "latin-1")
        anchor_points = self.compare_chroms(align_chrom, self.Method)
        if self.Method !="RI":
            anchor_points = self.check_rank(anchor_points, self.Ref_chrom)
            if len(anchor_points) ==0 or len(anchor_points) <5:
                return f"NO ENOUGH ALIGNMENT POINTS FOR {lst_element}"
            print(lst_element, "landmarks:", len(anchor_points))
        else:
            align_chrom = anchor_points
        lst_element = lst_element.split(".")[0]
        if self.Method == "DISCO":
            anchor_points.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.similarity}_{lst_element}_anchors.txt", sep = "\t")
        elif self.Method == "BiPACE2D":
            anchor_points.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.D1}_{self.D2}_{self.similarity}_{lst_element}_anchors.txt", sep = "\t")
        elif self.Method == "MSort":
            anchor_points.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.max_shift_1}_{self.max_shift_2}_{self.similarity}_{lst_element}_anchors.txt", sep = "\t")
        elif self.Method == "RI":
            pass
        elif self.Method == "TNTDA":
            anchor_points.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.w}_{self.distance}_{self.similarity}_{lst_element}_anchors.txt", sep = "\t")
        else:
            anchor_points.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.w}_{self.distance}_{lst_element}_anchors.txt", sep = "\t")
        if self.Method !="RI":
            align_chrom = self.shift_1stRT(anchor_points,self.Ref_chrom, align_chrom)
            align_chrom = self.shift_2ndRT(anchor_points,self.Ref_chrom, align_chrom)
        plt.scatter(self.Ref_chrom.loc[:,"1st Dimension Time (s)"],self.Ref_chrom.loc[:,"2nd Dimension Time (s)"], label = f"{self.Ref_chrom_name}_ref")
        plt.scatter(align_chrom.loc[:,"1st Dimension Time (s)"],align_chrom.loc[:,"2nd Dimension Time (s)"], label = f"{lst_element}_aligned")
        plt.legend(loc="upper left")
        if not os.path.isdir(self.Result_folder):
            os.makedirs(self.Result_folder)
        plt.draw()
        if self.Method == "DISCO":
            plt.savefig(f"{self.Result_folder}{os.sep}{self.Method}_{self.similarity}_{lst_element}.png", dpi=300)
            plt.close()
            align_chrom = align_chrom.set_index("Name")
            align_chrom.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.similarity}_{lst_element}.txt", sep = "\t")
        elif self.Method == "BiPACE2D":
            plt.savefig(f"{self.Result_folder}{os.sep}{self.Method}_{self.D1}_{self.D2}_{self.similarity}_{lst_element}.png",dpi=300)
            plt.close()
            align_chrom = align_chrom.set_index("Name")
            align_chrom.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.D1}_{self.D2}_{self.similarity}_{lst_element}.txt", sep = "\t")
        elif self.Method == "MSort":
            plt.savefig(f"{self.Result_folder}{os.sep}{self.Method}_{self.max_shift_1}_{self.max_shift_2}_{self.similarity}_{lst_element}.png",dpi=300)
            plt.close()
            align_chrom = align_chrom.set_index("Name")
            align_chrom.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.max_shift_1}_{self.max_shift_2}_{self.similarity}_{lst_element}.txt", sep = "\t")
        elif self.Method == "RI":
            plt.savefig(f"{self.Result_folder}{os.sep}{self.Method}_{lst_element}.png",dpi=300)
            plt.close()
            align_chrom = align_chrom.set_index("Name")
            align_chrom.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{lst_element}.txt", sep = "\t")
        elif self.Method == "TNTDA":
            plt.savefig(f"{self.Result_folder}{os.sep}{self.Method}_{self.w}_{self.distance}_{self.similarity}_{lst_element}.png", dpi=300)
            plt.close()
            align_chrom = align_chrom.set_index("Name")
            align_chrom.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.w}_{self.distance}_{self.similarity}_{lst_element}.txt", sep = "\t")
        else:
            plt.savefig(f"{self.Result_folder}{os.sep}{self.Method}_{self.w}_{lst_element}.png", dpi=300)
            plt.close()
            align_chrom = align_chrom.set_index("Name")
            align_chrom.to_csv(f"{self.Result_folder}{os.sep}{self.Method}_{self.w}_{lst_element}.txt", sep = "\t")
        end_time = datetime.now()
        if self.multipool == True:
            return f'{self.align_chrom_name} Duration: {end_time - start_time}'
        else:
            print(f'{self.align_chrom_name} Duration: {end_time - start_time}')
        
    def compare_chroms(self,align_chrom, method) ->pd.DataFrame:
        """
        Process the align_chrom and find the anchor peaks for the following DA.

        Parameters
        ----------
        align_chrom : pd.DataFrame
            Aligned chromatogram
        method : str
            DA method

        Returns
        -------
        result_table :  pd.DataFrame
            Final table of detected anchor peaks

        """
        check_table = align_chrom
        result_table = pd.DataFrame(index= self.Ref_chrom.index, columns = ["match", "1st Dimension Time (s)", "2nd Dimension Time (s)"])
        if method =="DISCO":
            result_table = self.DISCO(align_chrom, result_table, check_table)
        elif method == "BiPACE2D":
            result_table = self.BiPACE2D(align_chrom, result_table, check_table)
        elif method == "MSort":
            result_table = self.MSort(align_chrom, result_table, check_table)
        elif method== "RI":
            result_table = self.RI(align_chrom, result_table, check_table)
        elif method == "TNTDA":
            result_table = self.TNTDA(align_chrom, result_table, check_table)
        else:
            result_table = self.PAM(align_chrom, result_table, check_table)
        return result_table
                    

    def merge_peak_entries(self, filename:str)-> pd.DataFrame:
        """
        The function formats the column names into the proper strings, drops N/A rows, and removes "saturated( xxxx )" string from Area/Height columns
        The function iterates over all rows, identifies and merges all peaks of the same chemical origin.
        The function searches in the nearest surroundings considering the direction of an elution drift and compares the spectra of all peaks
        found in the investigated elution area.

        Parameters
        ----------
        filename : str
            Name of the file.

        Returns
        -------
        new_dataframe : pd.DataFrame
            Returns formated dataframe for further processing.

        """
        file = pd.read_csv(f"{self.Source_folder}{os.sep}{filename}", sep = "\t", encoding='latin-1').reset_index()
        try:
            file["2nd Dimension Time (s)"] = file["2nd Dimension Time (s)"].str.replace(",", ".").astype("float")
        except AttributeError:
            print("Data type for the transformation is already in a correct type.")
            pass
        try:
            file =file.rename(columns={"Spectra": "Spectrum"})
        except KeyError:
            print("Variable names are already correct .")
            pass
        ## this part of code reported as a bug causing crashes, will be fixed in future update
        # if sum(file.isna().sum())>0:
        #     file = file.dropna(how="all", axis=0)
        #     print("zero values dropped")
        if file["Area"].apply(isinstance,args = [str]).any()==True:
            file['Area'] = file['Area'].map(lambda x: x.lstrip('saturated( ').rstrip(' )')).astype("float")
            file['Height'] = file['Height'].map(lambda x: x.lstrip('saturated( ').rstrip(' )')).astype("float")
        print(filename)
        for row in file.index:
            df = pd.DataFrame(file.loc[row, :]).transpose()
            samplemasses_main = {}
            spectrum1 = df.at[row,"Spectrum"].strip().split(" ")
            for indic1 in spectrum1:
                key =int(float(indic1.split(":")[0]))
                value= float(indic1.split(":")[1])
                if key >=31:
                   if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                       value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                   samplemasses_main.update({key : value})
            lower_band_1 = int(df["1st Dimension Time (s)"])
            upper_band_1 = int(df["1st Dimension Time (s)"]+11)
            lower_band_2 = float(df["2nd Dimension Time (s)"])
            upper_band_2 = float(df["2nd Dimension Time (s)"]- 0.5)
            focus_range = file[file["1st Dimension Time (s)"].between(lower_band_1,upper_band_1)]
            focus_range = focus_range[focus_range["2nd Dimension Time (s)"].between(upper_band_2,lower_band_2)]
            if len(focus_range) >1:
                for sub_row in focus_range.index:
                    main_df = df
                    compare_df = pd.DataFrame(focus_range.loc[sub_row, :]).transpose()
                    if main_df.iat[0,0] == compare_df.iat[0,0]:
                        continue
                    else:
                        result = self.compare_spectra(samplemasses_main, compare_df.at[sub_row,"Spectrum"].strip().split(" "))
                        if result >= self.minimal_merge:
                            file.at[sub_row,"index"] = file.at[row,"index"]
                        pass
        new_dataframe = self.mask_groupby(file)
        new_dataframe = new_dataframe.drop("index", axis =1)
        new_dataframe = new_dataframe.set_index("Name")
        new_dataframe = new_dataframe.sort_values(by="1st Dimension Time (s)")
        print(f"{self.Folder_merged}{os.sep}{filename}")
        new_dataframe.to_csv(f"{self.Folder_merged}{os.sep}{filename}", sep = "\t")
        if self.multipool == True:
            return f"{filename}\nOriginal shape of Dataframe: {file.shape}||| New shape if Dataframe: {new_dataframe.shape}"
        else:
            print(f"{filename}\nOriginal shape of Dataframe: {file.shape}||| New shape if Dataframe: {new_dataframe.shape}")

    def compare_spectra(self, samplemasses_main:dict, spectrum2:list)->float:
        """
        Compares spectra of two peaks. Note that the spectra must be in basic ChromaTOF format mass:intensity separated by a blank space

        Parameters
        ----------
        spectrum1 : str
            mass:intensity format separated by a blank space
        spectrum2 : str
            mass:intensity format separated by a blank space
            
        Returns
        -------
        float
            The result of Pearson/cosine correlation coefficient in [%].

        """
        samplemasses_compare = dict()
        for indic2 in spectrum2:
            key =int(float(indic2.split(":")[0]))
            value= float(indic2.split(":")[1])
            if key >=31:
                if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                    value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                samplemasses_compare.update({key : value})
        spectrum_table_spectra = pd.DataFrame.from_dict([samplemasses_main, samplemasses_compare]).fillna(0)
        _1st_array = spectrum_table_spectra.iloc[0, :].to_numpy().astype("float")
        _2nd_array = spectrum_table_spectra.iloc[1, :].to_numpy().astype("float")
        if self.Similarity == "Pearson":
            pearson = stats.pearsonr(_1st_array,_2nd_array )
            result = round((pearson[0]*100),2)
        elif self.Similarity == "DOT":
            result = dot(_1st_array, _2nd_array)/(norm(_1st_array)*norm(_2nd_array))*100
        else:
            print("No valid method for spectral similarity. Force quit.")
            quit()
        return result
    
    def mask_groupby(self,df:pd.DataFrame)->pd.DataFrame:
        """
        Goes through the every row (peak) in a chromatogram and checks whenever there are two or more peaks of same chemical origin. If yes, 
        function merges them. Areas and Heights is sumed across all rows. Retention indexes are weighted and recalculated, mass spectra are inherited from the most
        intensive point of the merged peak (The most pure spectra are expected at the peaks´ climax)

        Parameters
        ----------
        df : pd.DataFrame
            Chromatogram with unmerged peaks

        Returns
        -------
        new_dataframe : pd.DataFrame
            Chromatogram with merged peaks.

        """
        lst_of_indexes = []
        new_dataframe = pd.DataFrame()
        for idx in df["index"]:
            if idx not in lst_of_indexes:
                lst_of_indexes.append(idx)
        for indx in lst_of_indexes:
            mask = df["index"].values == indx
            working_df = df[mask]
            if working_df.shape[0]==1:
                new_dataframe = pd.concat([new_dataframe, working_df], axis = 0, ignore_index=True, sort=False)
                continue
            else:
                repre_df = pd.DataFrame(columns = ["index","Name","1st Dimension Time (s)","2nd Dimension Time (s)", "Area", "Height", "Spectrum" ])
                working_df = working_df.fillna(0)
                area = working_df["Area"].sum()
                max_area = working_df["Area"].idxmax()
                height = working_df["Height"].sum()
                _1stRT = round((1/area)*sum(working_df["Area"]*working_df["1st Dimension Time (s)"]))
                _2ndRT = (1/area)*sum(working_df["Area"]*working_df["2nd Dimension Time (s)"])
                repre_df.loc [0] = [working_df.iat[0,0],working_df.iat[0,1],_1stRT,_2ndRT,area,height,working_df.at[max_area,"Spectrum" ]]
                new_dataframe = pd.concat([new_dataframe, repre_df], axis = 0, ignore_index=True, sort=False)
                pass
        del lst_of_indexes
        return new_dataframe
    
    def DISCO(self,align_chrom:pd.DataFrame, result_table:pd.DataFrame, check_table:pd.DataFrame)-> pd.DataFrame:
        """
        Based on Wang, B., et al., DISCO: Distance and Spectrum Correlation Optimization Alignment for Two-Dimensional Gas Chromatography Time-of-Flight Mass Spectrometry-Based Metabolomics. Analytical Chemistry, 2010. 82(12): p. 5069-5081.
        DOI: 10.1021/ac100064b
        To speed up processing speed, the algorithm considers only 20% of total rows with the least distance for each ref_table row.
        This parameter can be changed by bound = int((len(align_chrom)/100)*20) argument (line 440)

        Parameters
        ----------
        ref_chrom : pd.DataFrame
            Referential chromatogram
        align_chrom : pd.DataFrame
            Aligned chromatogram
        similarity : int
            Mass spectrum similarity treshold
        result_table : pd.DataFrame
            Template for the anchor peaks table
        check_table : pd.DataFrame
            Copy of original aligned chromatogram

        Returns
        -------
        result_table :  pd.DataFrame
            Final table of detected anchor peaks

        """
        ref_chrom = self.calculate_zscore_times(self.Ref_chrom)
        align_chrom = self.calculate_zscore_times(align_chrom)
        for i in ref_chrom.index:
            landmark_coords = (ref_chrom.at[i, "Z score 1st RT"],ref_chrom.at[i, "Z score 2nd RT"])
            align_chrom["Euclidian"] = np.sqrt(((landmark_coords[0]-align_chrom["Z score 1st RT"])**2)+((landmark_coords[1]-align_chrom["Z score 2nd RT"])**2))
            bound = int((len(align_chrom)/100)*20)
            sorted_df = (align_chrom.sort_values(by = "Euclidian", ascending = True)).iloc[0:bound]
            comp_results = {}
            for e in sorted_df.index:
                samplemasses_main = {}
                spectrum1 = ref_chrom.at[i,"Spectrum"].strip().split(" ")
                for indic1 in spectrum1:
                   key =int(float(indic1.split(":")[0]))
                   value= float(indic1.split(":")[1])
                   if key >=31:
                       if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                           value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                       samplemasses_main.update({key : value})
                   else:continue
                result = self.compare_spectra(samplemasses_main, sorted_df.at[e,"Spectrum"].strip().split(" "))
                if result>= self.similarity:
                    comp_results[e] = result
            try:
                max_value = max(comp_results, key = comp_results.get)
            except ValueError:
                max_value = "Not Found"
            if max_value != "Not Found":
                align_chrom = align_chrom.drop(max_value)
            result_table.at[i, "match"] = max_value
            if max_value != "Not Found":
                result_table.at[i, "1st Dimension Time (s)"] = check_table.at[max_value, "1st Dimension Time (s)"]
                result_table.at[i, "2nd Dimension Time (s)"] = check_table.at[max_value, "2nd Dimension Time (s)"]
            else:
                result_table.at[i, "match"] = "Not Found"
        result_table = result_table[result_table["match"] != "Not Found"]
        return result_table
    
    def BiPACE2D(self, align_chrom:pd.DataFrame, result_table:pd.DataFrame, check_table:pd.DataFrame)->pd.DataFrame:
        """
        The orginal algorithm:
        Hoffmann, N., et al., BiPACE 2D—graph-based multiple alignment for comprehensive 2D gas chromatography-mass spectrometry. Bioinformatics, 2013. 30(7): p. 988-995.
        DOI: https://doi.org/10.1093/bioinformatics/btt738
        This algorithm is the modification of the original one and was developed for the reaserch of human scent signature at UCT Prague.
    
        The algorithm takes 5 arguments: RT1 tolerance parameter, RT2 tolerance parameter, T1 critical value, T2 critical value and similarity critical value.
        Setting of T1, T2 influences the decision whether the algorithm should proceed with peak comparison or if the two comparised peaks are too far away.

        Parameters
        ----------
        align_chrom : pd.DataFrame
            Aligned chromatogram.
        result_table : pd.DataFrame
            Template for the anchor peaks table.
        check_table : pd.DataFrame
            Copy of original aligned chromatogram.
        D1 : float
            RT1 tolerance parameter.
        D2 : float
            RT2 tolerance parameter.
        T1 : float
            1st dimension treshold parameter.
        T2 : float
            1st dimension treshold parameter.

        Returns
        -------
        result_table : pd.DataFrame
            Final table of detected anchor peaks.

        """
        for i in self.Ref_chrom.index: 
          comp_results = {}
          for e in align_chrom.index:
              samplemasses_main = {}
              spectrum1 = self.Ref_chrom.at[i,"Spectrum"].strip().split(" ")
              for indic1 in spectrum1:
                 key =int(float(indic1.split(":")[0]))
                 value= float(indic1.split(":")[1])
                 if key >=31:
                     if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                         value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                     samplemasses_main.update({key : value})
                 else:continue
              rt1, rt2 = self.Ref_chrom.at[i,"1st Dimension Time (s)"], self.Ref_chrom.at[i,"2nd Dimension Time (s)"]
              comp_rt1, comp_rt2 = align_chrom.at[e,"1st Dimension Time (s)"], align_chrom.at[e,"2nd Dimension Time (s)"]
              if math.exp((-(rt1-comp_rt1)**2)/(2*self.D1)**2) >=self.T1 and math.exp((-(rt2-comp_rt2)**2)/(2*self.D2)**2) >= self.T2:
                  comp_spectra = self.compare_spectra(samplemasses_main, align_chrom.at[e,"Spectrum"].strip().split(" "))
                  if comp_spectra >= self.similarity:
                      result = comp_spectra*math.exp(((-(rt1-comp_rt1))**2/(2*self.D1)**2))*math.exp((-(rt2-comp_rt2))**2/(2*self.D2)**2)
                      comp_results[e] = result
                  else: continue
              else: continue
          try:
              max_value = max(comp_results, key = comp_results.get)
          except ValueError:
              max_value = "Not Found"
          if max_value != "Not Found":
              align_chrom = align_chrom.drop(max_value)
              result_table.at[i, "match"] = max_value
          if max_value != "Not Found":
              result_table.at[i, "1st Dimension Time (s)"] = check_table.at[max_value, "1st Dimension Time (s)"]
              result_table.at[i, "2nd Dimension Time (s)"] = check_table.at[max_value, "2nd Dimension Time (s)"]
          else:
              result_table.at[i, "match"] = "Not Found"
        result_table = result_table[result_table["match"] != "Not Found"]
        return result_table    
    
    def MSort(self,align_chrom:pd.DataFrame, result_table:pd.DataFrame, check_table:pd.DataFrame)->pd.DataFrame:
        """
        Oh, Cheolhwan, et al. "Comprehensive two-dimensional gas chromatography/time-of-flight mass spectrometry peak sorting algorithm." Journal of chromatography A 1179.2 (2008): 205-215.
        DOI: https://doi.org/10.1016/j.chroma.2007.11.101
        This algorithm is the modification of the original one and was developed for the reaserch of human scent signature at UCT Prague.
        
        Parameters max_shift_1 and max_shift_2 are the max allowed deviation of peaks in 1st and 2nd dimension (please note that the scale of 1st and 2nd differes)
        The similarity parameter is the threshold for the mass spectra similarity of compared peaks. 
        
        Parameters
        ----------
        align_chrom : pd.DataFrame
            Aligned chromatogram.
        result_table : pd.DataFrame
            Template for the anchor peaks table.
        check_table : pd.DataFrame
            Copy of original aligned chromatogram.
        Returns
        -------
        result_table : pd.DataFrame
            Final table of detected anchor peaks.

        """
        for i in self.Ref_chrom.index:
            samplemasses_main = {}
            spectrum1 = self.Ref_chrom.at[i,"Spectrum"].strip().split(" ")
            for indic1 in spectrum1:
               key =int(float(indic1.split(":")[0]))
               value= float(indic1.split(":")[1])
               if key >=31:
                   if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                       value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                   samplemasses_main.update({key : value})
               else:continue
            lower_band_1 = int(self.Ref_chrom.at[i,"1st Dimension Time (s)"]-self.max_shift_1)
            upper_band_1 = int(self.Ref_chrom.at[i,"1st Dimension Time (s)"]+self.max_shift_1)
            lower_band_2 = float(self.Ref_chrom.at[i,"2nd Dimension Time (s)"]-self.max_shift_2)
            upper_band_2 = float(self.Ref_chrom.at[i,"2nd Dimension Time (s)"]+self.max_shift_2)
            focus_range = align_chrom[align_chrom["1st Dimension Time (s)"].between(lower_band_1,upper_band_1)]
            focus_range = focus_range[focus_range["2nd Dimension Time (s)"].between(lower_band_2,upper_band_2)]
            if focus_range.shape[0] >0:
                comp_results = {}
                for e in focus_range.index:
                    result = self.compare_spectra(samplemasses_main, focus_range.at[e,"Spectrum"].strip().split(" "))
                    if result>= self.similarity:
                        comp_results[e] = result
                try:
                    max_value = max(comp_results, key = comp_results.get)
                except ValueError:
                    max_value = "Not Found"
                if max_value != "Not Found":
                    align_chrom = align_chrom.drop(max_value)
                result_table.at[i, "match"] = max_value
                if max_value != "Not Found":
                    result_table.at[i, "1st Dimension Time (s)"] = check_table.at[max_value, "1st Dimension Time (s)"]
                    result_table.at[i, "2nd Dimension Time (s)"] = check_table.at[max_value, "2nd Dimension Time (s)"]
            else:
                result_table.at[i, "match"] = "Not Found"
        result_table = result_table[result_table["match"] != "Not Found"]
        return result_table
  
    def PAM(self,align_chrom:pd.DataFrame, result_table:pd.DataFrame, check_table:pd.DataFrame)-> pd.DataFrame:
        """
        Based on Kim, Seongho, et al. "An optimal peak alignment for comprehensive two-dimensional gas chromatography mass spectrometry using mixture similarity measure." Bioinformatics 27.12 (2011): 1660-1666.
        DOI: https://doi.org/10.1093/bioinformatics/btr188
        This algorithm is the modification of the original one and was developed for the reaserch of human scent signature at UCT Prague.
        
        The w parameter is crucial for the calculation of mixture similarity: w/(1+Distance)*(1-w)*similarity
        
        Parameters
        ----------
        align_chrom : pd.DataFrame
            Aligned chromatogram.
        result_table : pd.DataFrame
            Template for the anchor peaks table.
        check_table : pd.DataFrame
            Copy of original aligned chromatogram.
        w : the weight for the mixture similarity measure
        
        Returns
        -------
        result_table :  pd.DataFrame
            Final table of detected anchor peaks
        """
        ref_chrom = self.Ref_chrom.reset_index()
        for row in ref_chrom.index:
            samplemasses_main = {}
            spectrum1 = self.Ref_chrom.at[row,"Spectrum"].strip().split(" ")
            for indic1 in spectrum1:
               key =int(float(indic1.split(":")[0]))
               value= float(indic1.split(":")[1])
               if key >=31:
                   if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                       value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                   samplemasses_main.update({key : value})
               else:continue
            t1, t2= ref_chrom.loc[row, "1st Dimension Time (s)"],ref_chrom.at[row, "1st Dimension Time (s)"]
            if self.distance == "Canberra":
                align_chrom["Dist"] = (abs(t1-align_chrom["1st Dimension Time (s)"])/abs(t1+align_chrom["1st Dimension Time (s)"]))+(abs(t2-align_chrom["2nd Dimension Time (s)"])/abs(t2+align_chrom["2nd Dimension Time (s)"])) ###Canberra
            elif self.distance == "Euklidian":
                align_chrom["Dist"] = np.sqrt((t1-align_chrom["1st Dimension Time (s)"])**2+(t2-align_chrom["2nd Dimension Time (s)"])**2) ###Euklid
            else:
                print("No valid method for calculating the distance. Force quit.")
                quit()
            align_chrom["Similarity"] = [self.compare_spectra(samplemasses_main, (align_chrom.at[x,"Spectrum"]).strip().split(" ")) for x in align_chrom.index]
            align_chrom["Result"] = (self.w/(1+align_chrom["Dist"]))+(1-self.w)*align_chrom["Similarity"]
            idx = align_chrom["Result"].idxmax()
            result_table.at[row, "match"] = idx
            result_table.at[row, "1st Dimension Time (s)"] = check_table.at[idx, "1st Dimension Time (s)"]
            result_table.at[row, "2nd Dimension Time (s)"] = check_table.at[idx, "2nd Dimension Time (s)"]
        return result_table
    
    def RI(self,align_chrom:pd.DataFrame, result_table:pd.DataFrame, check_table:pd.DataFrame)-> pd.DataFrame:
        """
        This algorithm was developed at UCT Prague for the purpose of the reaserch of human scent signature.
        The algorithm requires an additional input described in software documentation.
        The fixed anchor points are identified by the user. The retention times are recalculated to retention indexes (note the name of columns
                                                                                                                       remain unchanged)

         Parameters
         ----------
         align_chrom : pd.DataFrame
             Aligned chromatogram.
         result_table : pd.DataFrame
             Template for the anchor peaks table.
         check_table : pd.DataFrame
             Copy of original aligned chromatogram.
         w : the weight for the mixture similarity measure
         
         Returns
         -------
         result_table :  pd.DataFrame
             Final table of detected anchor peaks

        """
        table = self.align_chrom_name.split("_")[-2]+"_"+self.align_chrom_name.split("_")[-1]
        transformed_df_1 = pd.DataFrame(columns = align_chrom.columns)
        transformed_df_2= pd.DataFrame(columns = align_chrom.columns)
        ref_df= pd.read_csv(f"metadata{os.sep}ref_data{os.sep}{table}", sep="\t", header=0, index_col=0)
        lst_1stD = [int(ref_df.loc[x,"R.T. (s)"].split(",")[0]) for x in ref_df.index if "Ref1" in x]
        lst_1stD = sorted(lst_1stD)
        lst_2ndD = [float(ref_df.loc[x,"R.T. (s)"].split(",")[1]) for x in ref_df.index if "Ref2" in x or "Ref1_8" in x]
        lst_2ndD = sorted(lst_2ndD)
        if len(lst_1stD) != 13 or len(lst_2ndD) != 5:
            print("Error!", self.align_chrom_name)
        df_low_1stD = align_chrom[align_chrom['1st Dimension Time (s)'].values<lst_1stD[0]]
        if len(df_low_1stD) >0:
            df_low_1stD.loc[:,'1st Dimension Time (s)']= 10+(np.log10(df_low_1stD.loc[:,'1st Dimension Time (s)'])-np.log10(lst_1stD[0]))/(np.log10(lst_1stD[1])-np.log10(lst_1stD[0]))
            transformed_df_1 = transformed_df_1.append(df_low_1stD, ignore_index = False)
            del df_low_1stD
        for loop in range(1,len(lst_1stD)):
            low = lst_1stD[loop-1]
            high = lst_1stD[loop] 
            df_foc = align_chrom[align_chrom['1st Dimension Time (s)'].between(low+0.00001,high)]
            df_foc.loc[:,'1st Dimension Time (s)'] = (10+loop)+(np.log10(df_foc.loc[:,'1st Dimension Time (s)'])-np.log10(low))/(np.log10(high)-np.log10(low))
            transformed_df_1 = transformed_df_1.append(df_foc, ignore_index = False)
            del df_foc
        df_high_1stD = align_chrom[align_chrom['1st Dimension Time (s)']>lst_1stD[-1]]
        if len(df_high_1stD) >0:
            df_high_1stD.loc[:,'1st Dimension Time (s)']= 21+(np.log10(df_high_1stD.loc[:,'1st Dimension Time (s)'])-np.log10(lst_1stD[-2]))/(np.log10(lst_1stD[-1])-np.log10(lst_1stD[-2]))
            transformed_df_1 = transformed_df_1.append(df_high_1stD, ignore_index = False)
            del df_high_1stD 
        df_low_2ndD = transformed_df_1[transformed_df_1['2nd Dimension Time (s)'].values<lst_2ndD[0]]
        if len(df_low_2ndD) >0:
            df_low_2ndD.loc[:,'2nd Dimension Time (s)']= 1+(np.log10(df_low_2ndD.loc[:,'2nd Dimension Time (s)'])-np.log10(lst_2ndD[0]))/(np.log10(lst_2ndD[1])-np.log10(lst_2ndD[0]))
            transformed_df_2 = transformed_df_2.append(df_low_2ndD, ignore_index = False)
            del df_low_2ndD
        for loop in range(0,len(lst_2ndD)-1):
            low = lst_2ndD[loop]
            high = lst_2ndD[loop+1] 
            df_foc = transformed_df_1[transformed_df_1['2nd Dimension Time (s)'].between(low+0.00001,high)]
            df_foc.loc[:,'2nd Dimension Time (s)'] = (1+loop)+(np.log10(df_foc.loc[:,'2nd Dimension Time (s)'])-np.log10(low))/(np.log10(high)-np.log10(low))
            transformed_df_2 = transformed_df_2.append(df_foc, ignore_index = False)
            del df_foc
        df_high_2ndD = transformed_df_1[transformed_df_1['2nd Dimension Time (s)']>lst_2ndD[-1]]
        if len(df_high_2ndD) >0:
            df_high_2ndD.loc[:,'2nd Dimension Time (s)']= 4+(np.log10(df_high_2ndD.loc[:,'2nd Dimension Time (s)'])-np.log10(lst_2ndD[-2]))/(np.log10(lst_2ndD[-1])-np.log10(lst_2ndD[-2]))
            transformed_df_2 = transformed_df_2.append(df_high_2ndD, ignore_index = False)
            del df_high_2ndD
        return transformed_df_2
    
    def TNTDA(self,align_chrom:pd.DataFrame, result_table:pd.DataFrame, check_table:pd.DataFrame)-> pd.DataFrame:
        """
        Derived algorithm which combines discrimination parameters of DISCO and selection rule of PAM algorithm.

        Parameters
        ----------
        align_chrom : pd.DataFrame
            Aligned chromatogram.
        result_table : pd.DataFrame
            Template for the anchor peaks table.
        check_table : pd.DataFrame
            Copy of original aligned chromatogram.
        w : the weight for the mixture similarity measure; w/(1+Distance)*(1-w)*similarity
        The similarity parameter is the threshold for the mass spectra similarity of compared peaks. 
        Canberra or Euclidian are available for calculating the distance between the datapoints. 
        
        
        Returns
        -------
        result_table : TYPE
            DESCRIPTION.

        """
        ref_chrom = self.calculate_zscore_times(self.Ref_chrom)
        align_chrom = self.calculate_zscore_times(align_chrom)
        for i in ref_chrom.index:
            t1, t2 = (ref_chrom.at[i, "Z score 1st RT"],ref_chrom.at[i, "Z score 2nd RT"])
            if self.distance == "Canberra":
                align_chrom["Dist"] = (abs(t1-align_chrom["Z score 1st RT"])/abs(t1+align_chrom["Z score 1st RT"]))+(abs(t2-align_chrom["Z score 2nd RT"])/abs(t2+align_chrom["Z score 2nd RT"])) ###Canberra
            elif self.distance == "Euklidian":
                align_chrom["Dist"] = np.sqrt((t1-align_chrom["Z score 1st RT"])**2+(t2-align_chrom["Z score 2nd RT"])**2) ###Euklid
            bound = int((len(align_chrom)/100)*20)
            sorted_df = (align_chrom.sort_values(by = "Dist", ascending = True)).iloc[0:bound]
            comp_results = {}
            for e in sorted_df.index:
                samplemasses_main = {}
                spectrum1 = ref_chrom.at[i,"Spectrum"].strip().split(" ")
                for indic1 in spectrum1:
                   key =int(float(indic1.split(":")[0]))
                   value= float(indic1.split(":")[1])
                   if key >=31:
                       if self.Transformation[0] != 1 or self.Transformation[1] != 0:
                           value = round((key**self.Transformation[1])*(value**self.Transformation[0]),2) 
                       samplemasses_main.update({key : value})
                   else:continue
                spectral_sim = self.compare_spectra(samplemasses_main, sorted_df.at[e,"Spectrum"].strip().split(" "))
                if spectral_sim>= self.similarity:
                    result = (self.w/(1+align_chrom.at[e,"Dist"]))+(1-self.w)*spectral_sim
                    comp_results[e] = result
            try:
                max_value = max(comp_results, key = comp_results.get)
            except ValueError:
                max_value = "Not Found"
            if max_value != "Not Found":
                align_chrom = align_chrom.drop(max_value)
            result_table.at[i, "match"] = max_value
            if max_value != "Not Found":
                result_table.at[i, "1st Dimension Time (s)"] = check_table.at[max_value, "1st Dimension Time (s)"]
                result_table.at[i, "2nd Dimension Time (s)"] = check_table.at[max_value, "2nd Dimension Time (s)"]
            else:
                result_table.at[i, "match"] = "Not Found"
        result_table = result_table[result_table["match"] != "Not Found"]
        return result_table
    
    def calculate_zscore_times(self,dataframe:pd.DataFrame)-> pd.DataFrame:
        """
        Recalculates original retention times to Z-scores. The result daraframe has n+2 columns
        "Z score 1st RT" and "Z score 2nd RT" columns are added at the end of the table (last two indexes)

        Parameters
        ----------
        dataframe : pd.DataFrame

        Returns
        -------
        dataframe : pd.DataFrame

        """
        if "Z score 1st RT" in dataframe.columns:
            pass
        else:
            dataframe.insert(len(dataframe.columns),"Z score 1st RT",0, True)
            dataframe.insert(len(dataframe.columns),"Z score 2nd RT",0, True)
            pass
        dataframe["Z score 1st RT"] = (dataframe["1st Dimension Time (s)"]-(dataframe["1st Dimension Time (s)"].mean()))/(dataframe["1st Dimension Time (s)"].std())
        dataframe["Z score 2nd RT"] = (dataframe["2nd Dimension Time (s)"]-(dataframe["2nd Dimension Time (s)"].mean()))/(dataframe["2nd Dimension Time (s)"].std())
        return dataframe
    
    def check_rank (self, anchor_points:pd.DataFrame, ref_table:pd.DataFrame)-> pd.DataFrame:
        """
        Checks the elution order of anchor points in both dimension via absolute retention time.

        Parameters
        ----------
        anchor_points : pd.DataFrame
            Dataframe of anchor point found in aligned chromatogram
        ref_table : pd.DataFrame
            Dataframe of anchor point found in referential chromatogram.

        Returns
        -------
        anchor_points : pd.DataFrame
            Updated dataframe of anchor point found in aligned chromatogram

        """
        anchor_points["1st RT_ref"] = [ref_table.at[x,"1st Dimension Time (s)"] for x in anchor_points.index]
        anchor_points["2nd RT_ref"] = [ref_table.at[x,"2nd Dimension Time (s)"] for x in anchor_points.index]
        anchor_points["Absolute ref RT"] = anchor_points["1st RT_ref"]+anchor_points["2nd RT_ref"]
        anchor_points["Absolute aligned RT"] = anchor_points["1st Dimension Time (s)"]+anchor_points["2nd Dimension Time (s)"]
        anchor_points.drop_duplicates(subset = ["Absolute ref RT","Absolute aligned RT"], keep="first", inplace=True)
        check_table = anchor_points
        anchor_points = anchor_points.sort_values(by = "Absolute ref RT", ascending = True)
        anchor_points["Rank_ref"] = range(0,len(anchor_points))
        anchor_points = anchor_points.sort_values(by = "Absolute aligned RT", ascending = True)
        anchor_points["Rank_match"] = range(0,len(anchor_points))
        anchor_points["delta"] = abs( anchor_points["Rank_ref"]-anchor_points["Rank_match"])
        check = (anchor_points["delta"] == 0).all()
        while check == False:
            max_val = anchor_points["delta"].max()
            mask = anchor_points["delta"].values == max_val
            working_df = anchor_points[mask]
            drop_rows = [x for x in working_df.index]
            anchor_points = anchor_points.drop(drop_rows[0])
            check_table = check_table.drop(drop_rows[0])
            anchor_points = check_table
            anchor_points = anchor_points.sort_values(by = "Absolute ref RT", ascending = True)
            anchor_points["Rank_ref"] = range(0,len(anchor_points))
            anchor_points = anchor_points.sort_values(by = "Absolute aligned RT", ascending = True)
            anchor_points["Rank_match"] = range(0,len(anchor_points))
            anchor_points["delta"] = abs( anchor_points["Rank_ref"]-anchor_points["Rank_match"])
            check = (anchor_points["delta"] == 0).all()
        anchor_points.drop(["Rank_ref","Rank_match", "delta"], axis = 1, inplace = True)
        anchor_points = anchor_points.sort_values(by ="1st Dimension Time (s)", ascending = True)
        return anchor_points
    
    def shift_1stRT(self,anchor_points:pd.DataFrame,ref_chrom:pd.DataFrame, align_chrom:pd.DataFrame)->pd.DataFrame:
        """
        Function takes the list of anchor points and corrects the 1st dimension retention time shift.

        Parameters
        ----------
        anchor_points : pd.DataFrame
            Anchor points with ID links to referential and aligned chromatogram.
        ref_chrom : pd.DataFrame
            Referential chromatogram.
        align_chrom : pd.DataFrame
            Aligned chromatogram.

        Returns
        -------
        transformed_df : pd.DataFrame
            Dataframe with corrected 1st dimensional shift

        """
        anchor_points.drop_duplicates(subset = ["1st Dimension Time (s)","1st RT_ref"], keep="first", inplace=True)
        transformed_df = pd.DataFrame(columns = align_chrom.columns)
        ref_points = [ref_chrom.at[x,"1st Dimension Time (s)"] for x in anchor_points.index]
        ref_points.sort()
        breaking_points = [align_chrom.at[x,"1st Dimension Time (s)"] for x in anchor_points["match"].values]
        breaking_points.sort()
        crit_values_target = breaking_points[0], breaking_points[-1]
        database_foc_low = align_chrom[align_chrom["1st Dimension Time (s)"].values <=crit_values_target[0]]
        if len(database_foc_low)>0:
            ref_value = ref_points[0]
            database_foc_low.loc[:,"1st Dimension Time (s)"] = (database_foc_low.loc[:,"1st Dimension Time (s)"]/crit_values_target[0])*ref_value
            transformed_df = pd.concat([transformed_df,database_foc_low], ignore_index = False)
            del database_foc_low
        database_foc_high = align_chrom[align_chrom["1st Dimension Time (s)"].values >crit_values_target[1]]
        if len(database_foc_high)>0:
            ref_value = ref_points[-1]
            database_foc_high.loc[:,"1st Dimension Time (s)"] = (database_foc_high["1st Dimension Time (s)"]/crit_values_target[1])*ref_value
            transformed_df = pd.concat([transformed_df,database_foc_high], ignore_index = False)
            del database_foc_high
        loops = len(breaking_points)-1
        for loop in range(loops):
            low_value_target = breaking_points[loop]
            high_value_target = breaking_points[loop+1]
            low_value_ref = ref_points[loop]
            high_value_ref = ref_points[loop+1]
            database_foc = align_chrom[align_chrom["1st Dimension Time (s)"].between(low_value_target+0.0001, high_value_target)]
            database_foc.loc[:,"1st Dimension Time (s)"] = ((database_foc["1st Dimension Time (s)"]-low_value_target)/(high_value_target-low_value_target))*(high_value_ref-low_value_ref)+low_value_ref
            transformed_df = pd.concat([transformed_df,database_foc], ignore_index = False)
            del database_foc
        transformed_df = transformed_df.sort_values(by ="1st Dimension Time (s)" )
        del align_chrom
        return transformed_df

    def shift_2ndRT(self,anchor_points:pd.DataFrame,ref_chrom:pd.DataFrame, align_chrom:pd.DataFrame)->pd.DataFrame:
        """
        Function takes the list of anchor points and corrects the 2nd dimension retention time shift.

        Parameters
        ----------
        anchor_points : pd.DataFrame
            Anchor points with ID links to referential and aligned chromatogram.
        ref_chrom : pd.DataFrame
            Referential chromatogram.
        align_chrom : pd.DataFrame
            Aligned chromatogram.

        Returns
        -------
        transformed_df : pd.DataFrame
            Dataframe with corrected 2nd dimensional shift

        """
        anchor_points.drop_duplicates(subset = ["2nd Dimension Time (s)","2nd RT_ref"], keep="first", inplace=True)
        anchor_points = anchor_points.sort_values(by ="2nd Dimension Time (s)", ascending = True)
        transformed_df = pd.DataFrame(columns = align_chrom.columns)
        ref_points = [ref_chrom.at[x,"2nd Dimension Time (s)"] for x in anchor_points.index]
        ref_points.sort()
        breaking_points = [align_chrom.at[x,"2nd Dimension Time (s)"] for x in anchor_points["match"].values]
        breaking_points.sort()
        crit_values_target = breaking_points[0], breaking_points[-1]
        database_foc_low = align_chrom[align_chrom["2nd Dimension Time (s)"].values <=crit_values_target[0]]
        if len(database_foc_low)>0:
            ref_value = ref_points[0]
            database_foc_low.loc[:,"2nd Dimension Time (s)"] = (database_foc_low["2nd Dimension Time (s)"]/crit_values_target[0])*ref_value
            transformed_df = transformed_df.append(database_foc_low, ignore_index = False)
            del database_foc_low
        database_foc_high = align_chrom[align_chrom["2nd Dimension Time (s)"].values >crit_values_target[1]]
        if len(database_foc_high)>0:
            ref_value = ref_points[-1]
            database_foc_high.loc[:,"2nd Dimension Time (s)"] = (database_foc_high["2nd Dimension Time (s)"]/crit_values_target[1])*ref_value
            transformed_df = transformed_df.append(database_foc_high, ignore_index = False)
            del database_foc_high
        loops = len(breaking_points)-1
        for loop in range(loops):
            low_value_target = breaking_points[loop]
            high_value_target = breaking_points[loop+1]
            low_value_ref = ref_points[loop]
            high_value_ref = ref_points[loop+1]
            if low_value_ref == high_value_ref:
                high_value_ref = ref_points[loop+1]+(high_value_target-low_value_target)
            database_foc = align_chrom[align_chrom["2nd Dimension Time (s)"].between(low_value_target+0.00001, high_value_target)]
            database_foc.loc[:,"2nd Dimension Time (s)"] = ((database_foc["2nd Dimension Time (s)"]-low_value_target)/(high_value_target-low_value_target))*(high_value_ref-low_value_ref)+low_value_ref
            transformed_df = transformed_df.append(database_foc, ignore_index = False)
            del database_foc
        transformed_df = transformed_df.sort_values(by ="1st Dimension Time (s)" )
        del align_chrom
        return transformed_df
if __name__ == '__main__':
    x = data_align("230815_LR_07_EUCS_Fr1", True, method = "DISCO", similarity="DOT", transform =(0.53,1.3), multipool=True).run_process()