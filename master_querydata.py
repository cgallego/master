# -*- coding: utf-8 -*-
"""
Master python script to run each module in sequence
Created on Tue May 20 11:17:08 2014


Arguments:
============
sys.argv[1] = input text file with one case per line in the minimum following format:
BenignNMaligNAnt	StudyNumber	DicomExamNumber	LesionID	StudyDate	SeriesID	BreastSide	PathReportID	PathoBreastSide	


@ author (C) Cristina Gallego, University of Toronto
--------------------------------------------------------------------
 """

import os, os.path
import sys
import string

from sys import argv, stderr, exit
import shlex, subprocess
import re

import numpy as np
import dicom
import psycopg2
import sqlalchemy as al
import sqlalchemy.orm

import pandas as pd
from query_database import *

from dictionaries import my_aet, hostID, local_port, clinical_aet, clinical_IP, clinical_port, remote_aet, remote_IP, remote_port
import dcmtk_routines as dcmtk

from inputs_init import *
from display import *
from segment import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
import pylab    
 
# convertion packages
import pandas.rpy.common as com
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as R
from rpy2.robjects import globalenv
from classifyCascade import *


def getScans(path_rootFolder, fileline, PatientID, StudyID, AccessionN, oldExamID):
    """
    run : getScans(path_rootFolder, PatientID, StudyID, AccessionN):

    Inputs
    ======
    path_rootFolder: (string)   Automatically generated based on the location of file 
    
    PatientID : (int)    MRN
    
    StudyID : (int)    CAD StudyID
    
    AccessionN : (int)  CAD AccessionN
    
    database : (bool) [True]   whether to check database for info about study.
    
    Output
    ====== 
    """
    try:
        dcmtk.check_MRI_MARTEL(data_loc, remote_aet, remote_port, remote_IP, local_port, PatientID, StudyID, AccessionN)
        if(oldExamID==False):
            dcmtk.pull_MRI_MARTEL(path_rootFolder, data_loc, remote_aet, remote_port, remote_IP, local_port, PatientID, StudyID, AccessionN, countImages=False)
        else:
            ExamID = fileline[4]
            dcmtk.pull_MRI_MARTELold(path_rootFolder, data_loc, remote_aet, remote_port, remote_IP, local_port, PatientID, StudyID, AccessionN, ExamID, countImages=False)
            
    except (KeyboardInterrupt, SystemExit):
        dcmtk.check_pacs(path_rootFolder, data_loc,  clinical_aet , clinical_port, clinical_IP, local_port, PatientID, StudyID, AccessionN)
        dcmtk.pull_pacs(path_rootFolder, data_loc, clinical_aet, clinical_port, clinical_IP, local_port, PatientID, StudyID, AccessionN)
    except (KeyboardInterrupt, SystemExit):
        print 'Unable to find study in MRI_MARTEL or AS0SUNB --- Abort'
        sys.exit()
        
    return
    
 
if __name__ == '__main__':
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
   
    # Open filename list
    file_ids = open(sys.argv[1],"r")
    init_flag=True
    
    for fileline in file_ids:
        # Get the line: StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID    image_pos_pat    image_ori_pat
        fileline = fileline.split()
        cond = fileline[0] 
        StudyID = fileline[1]  
        DicomExamNumber = fileline[2]
        Lesions_id = fileline[3]
        dateID = fileline[4]
        SeriesID = fileline[5] # corresponds to dynamic sequence;
        
        #############################
        ###### 1) Retrieving Images from Research PACS
        #############################
        print "Retrieving Scans to local drive..."
        #getScans(path_rootFolder, fileline, PatientID, StudyID, AccessionN, oldExamID=False)
        
        #############################
        ###### 2) Querying Research database for clinical, pathology, radiology data
        #############################
        print "Executing SQL connection..."
        # Format query StudyID
        if (len(StudyID) >= 4 ): fStudyID=StudyID
        if (len(StudyID) == 3 ): fStudyID='0'+StudyID
        if (len(StudyID) == 2 ): fStudyID='00'+StudyID
        if (len(StudyID) == 1 ): fStudyID='000'+StudyID
           
        # Format query redateID
        redateID = dateID[0:4]+'-'+dateID[4:6]+'-'+dateID[6:8]
    
        # perform query        
        queryData = Query()
        queryData.queryDatabase(fStudyID, redateID)       
        rowCase=["0", "0"]
        rowCase = int(raw_input('pick row (0-n): '))
        
        # recollect pathologies
        queryData.d1['is_insitu'] = pd.Series(True, index=queryData.d1)
        queryData.d1['is_invasive'] = pd.Series(True, index=queryData.d1)
        queryData.d1['Diagnosis'] = pd.Series(True, index=queryData.d1)
        queryData.d1['BenignNMaligNAnt'] = pd.Series(True, index=queryData.d1)
        queryData.d1['labels'] = pd.Series(True, index=queryData.d1)
        ansLesion = array((raw_input('Enter: is_insitu?: is_invasive?: ')).split()).astype(bool)

        #slice data, get only 1 record        
        dataCase = pd.Series( queryData.d1.loc[rowCase,:] )
        dataCase['is_insitu'] =  ansLesion[0]
        dataCase['is_invasive'] =  ansLesion[1]
        
        ansDiag=["diagnosis"]
        ansDiag = str(raw_input('Dignosis: '))
        dataCase['Diagnosis'] =  ansDiag
        dataCase['BenignNMaligNAnt'] =  cond[:-1]
        dataCase['labels'] =  cond
        
        if(init_flag): 
            casesFrame = pd.DataFrame(columns=queryData.d1.columns)
            init_flag=False
            
        #############################
        ###### Finish tidying up and save to file
        ## append collection of cases
        #############################
        casesFrame = casesFrame.append(dataCase) # 20
        casesFrame['id']=fStudyID
        casesFrame.set_index('id',inplace=False)
             
    file_ids.close()
    # end of line
    casesFrame.to_csv('casesFrames_queried.csv')   