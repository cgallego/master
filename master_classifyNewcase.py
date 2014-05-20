# -*- coding: utf-8 -*-
"""
Master python script to run each module in sequence

Arguments:
============
sys.argv[1] = input text file with one case per line in the following format:
BenignNMaligNAnt	StudyNumber	DicomExamNumber	LesionID	StudyDate	SeriesID	BreastSide	PathReportID	PathoBreastSide	

Created on Tue Apr 08 14:20:35 2014
@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
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
    init_flag=1
    
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
        ###### 3) Extractfeatures
        #############################
        ###### Start by Loading 
        print "Start by loading volumes..."
        load = Inputs_init()
        data_loc='Z:\Cristina\MassNonmass'+os.sep+cond[:-1]
        
        [series_path, phases_series, lesionID_path] = load.readVolumes(data_loc, StudyID, DicomExamNumber, SeriesID, Lesions_id)
        print "Path to series location: %s" % series_path 
        print "List of pre and post contrast volume names: %s" % phases_series
        print "Path to lesion segmentation: %s" % lesionID_path
        
        print "\n Load Segmentation..."
        lesion3D = load.loadSegmentation(lesionID_path)
        print "Data Structure: %s" % lesion3D.GetClassName()
        print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
        print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())
        
        print "\n Visualize volumes..."
        # Create only 1 display
        loadDisplay = Display()
        lesion3D_mesh = loadDisplay.addSegment(lesion3D, (0,1,0), interact=False)
        loadDisplay.visualize(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, sub=True, postS=4, interact=True)
        
        # Get z slice
        LesionZslice = loadDisplay.zImagePlaneWidget.GetSliceIndex()

        #############################
        # 4) Create Segmentation of lesion. Comment out if not needed ( define seededlesion3D = lesion3D  )
        #############################
        createSegment = Segment()
        print "\n Displaying picker for lesion segmentation"
        seeds = loadDisplay.display_pick(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, 4, LesionZslice)
        
        seededlesion3D = createSegment.segmentFromSeeds(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, seeds, loadDisplay.iren1, loadDisplay.xImagePlaneWidget, loadDisplay.yImagePlaneWidget,  loadDisplay.zImagePlaneWidget)
        seededlesion3D_mesh = loadDisplay.addSegment(seededlesion3D, (0,0,1), interact=True)
        loadDisplay.picker.RemoveAllObservers()

        # save it to file	             
        createSegment.saveSegmentation(lesionID_path, seededlesion3D)                
        
        #############################
        ###### Extract Dynamic features
        #############################
        print "\n Extract Dynamic contour features..."
        loadDynamic = Dynamic()
        dynamicfeatures_contour = loadDynamic.extractfeatures_contour(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, seededlesion3D)
        print "\n=========================================="
        print dynamicfeatures_contour
        
        print "\n Extract Dynamic inside features..."
        dynamicfeatures_inside = loadDynamic.extractfeatures_inside(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, seededlesion3D)
        print dynamicfeatures_inside
        print "\n=========================================="
        
        #############################
        ###### Extract Morphology features
        #############################
        print "\n Extract Morphology features..."
        loadMorphology = Morphology()
        morphofeatures = loadMorphology.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, seededlesion3D)
        print "\n=========================================="
        print morphofeatures
        print "\n=========================================="

        #############################        
        ###### Extract Texture features
        #############################
        print "\n Extract Texture features..."
        loadTexture = Texture()
        texturefeatures = loadTexture.extractfeatures(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, series_path, phases_series, seededlesion3D, loadMorphology.VOI_efect_diameter, loadMorphology.lesion_centroid )
        print "\n=========================================="
        print texturefeatures
        print "\n=========================================="
        
        # deal with closing windows, plots, renders, actors
        pylab.close('all')
        loadDisplay.renderer1.RemoveActor(loadDisplay.actor_mesh)
        loadDisplay.iren1.TerminateApp()
        loadDisplay.renWin1.Finalize()
        
        #############################
        ###### Finish tidying up and save to file
        ## append collection of cases
        #############################
        casesFrame = casesFrame.append(dataCase) # 20
        casesFrame['id']=fStudyID
        casesFrame.set_index('id',inplace=False)
        
        dynamicfeatures_contour['id']=fStudyID
        dynamicfeatures_contour.set_index('id',inplace=False)
        
        casesFrame = pd.merge(casesFrame, dynamicfeatures_contour, on='id', how='inner')
        
        dynamicfeatures_inside['id']=fStudyID
        dynamicfeatures_inside.set_index('id',inplace=False)
        casesFrame = pd.merge(casesFrame, dynamicfeatures_inside, on='id', how='inner')
        
        morphofeatures['id']=fStudyID
        morphofeatures.set_index('id',inplace=False)
        casesFrame = pd.merge(casesFrame, morphofeatures, on='id', how='inner')
        
        texturefeatures['id']=fStudyID
        texturefeatures.set_index('id',inplace=False)
        casesFrame = pd.merge(casesFrame, texturefeatures, on='id', how='inner')
        
        # end of line
        os.chdir("Z:\Cristina\MassNonmass\codeProject\codeBase\extractFeatures\casesDatabase")
        casesFrame.to_csv('casesFrames_toclasify.csv')   
        
        #############################
        ## Classification stage: send features to classifier to generate new case prediction
        ## Cascade current classifier
        #############################
        classifier = classifyCascade()
        classifier.case_classifyCascade()
             
    file_ids.close()
    