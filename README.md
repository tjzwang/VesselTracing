# VesselTracing
  
Protocol:   
1. Calculate the intensity of objects or pixels along the length of blood vessels  
       Files:  
         intensity_distance_analysis_NEW.m (Use this file if data is correctly oriented. You can check orientation of data by looking at figures that appear. In figure 4, the red vessel should be visible inside the blue tau data) 
       
     intensity_distance_analysis_NEW_FLIPPED.m (Use if vessel is not inside tau data. This script flips data to reorient the vessel correctly. Check figure 4 to ensure flipping fixed the misalignment.)  
    
2. Bin data along and away from vessel:  
      Files:  
         For Intensity: Vessel_stat_folder_intensity_binning.m  
         For Objects: Vessel_stat_folder_object_binning.m  
  
NOTE: The folder "Necessary Functions" must be added to the path for scripts to run!  
  
MATLAB Code Descriptions:  
  
1.    intensity_distance_analysis_NEW.m    
      This script calculates the distance of pixels along and away from blood vessels.  
      Inputs:  
          Data:  
             1. Tau Coordinates (Obtain this by exporting coordinates of original tau channel in FIJI. Should be a .txt file with XYZ and intensity values for each pixel).  
             2. Vessel Data (Coordinates for each blood vessel. Obtain these by loading Imaris vessel masks into FIJI and exporting coordinates for each vessel. Should be .txt files.)  
             3. Distance Transforms (Distance transform for each vessel. Obtain these first in Imaris, then refine with FIJI MACRO to reduce file size. Should be 1 per vessel. Ensure files are named in similar form to the vessel data so they are ordered correctly.)  
          Values:  
             1. Imaris Image Dimensions in Microns (Includes lower and upper bounds for X, Y, and Z.)  
             2. Fiji Image Dimensions in Pixel Count (Include X dimension, Y dimension, and # of Z-stacks for each input: tau data, vessel data, and distance transform data)  
             3. Export path, name, and file type for exported data files. (.txt is recommended)  
      Output Format:  
       Column: [  1    2    3         4                  5                      6           ]  
               [  X    Y    Z     Intensity     DistanceAlongVessel    DistanceFromSurface  ]
       Data exports as .txt or .xlsx (.txt is recommended)  
               
2.   intensity_distance_analysis_NEW_FLIPPED.m  
     Script completes same function asintensity_distance_analysis_NEW.m  
     Inputs:  
         Same inputs as intensity_distance_analysis_NEW.m. Only use if data is not correctly aligned in intensity_distance_analysis_NEW.m.  
     Output:  
         Same format as intensity_distance_analysis_NEW.m  
           
3.   Vessel_stat_folder_intensity_binning.m   
     Bins intensity data in groups based on distance along the vessel and distance away from the vessel.   
     Inputs:   
          1. Folder containing files with calculated pixel distances along and away from each blood vessel (Outputs from intensity_distance_analysis_NEW.m)  
     Values:   
          1. Max Distance From Vessel (in microns)  
          2. Min Distance From Vessel (in microns)  
          3. Distance from Vessel Resolution (in microns)  
          4. Binning Along Vessel Frame Width (in microns)  
     Output Format:   
      Column: [    1        2            3                4                5-(Max Distance From Vessel - Min Distance From Vessel  ]  
              [ Vessel#  Frame#  MinFrameDistance    MaxFrameDistance          Mean Intensity For Each Distance From Vessel        ]  
      Data exports as a .csv
            
4. Vessel_stat_folder_object_binning.m  
      Bins object data in groups based on distance along the vessel and distance away from the vessel.   
      Inputs: 
          1. Folder containing files with calculated object distances along and away from each blood vessel (Outputs from intensity_distance_analysis_NEW.m)  
     Values:   
          1. Max Distance From Vessel (in microns)  
          2. Min Distance From Vessel (in microns)  
          3. Distance from Vessel Resolution (in microns)   
          4. Binning Along Vessel Frame Width (in microns)   
     Output Format:   
      Column: [    1        2            3                4                5-(Max Distance From Vessel - Min Distance From Vessel  ]  
              [ Vessel#  Frame#  MinFrameDistance    MaxFrameDistance           (Density For Each Distance From Vessel)            ]  
      Data exports as a .csv  
   
5. Necessary Functions  
     This folder contains functions necessary for the other scripts to run. These functions must be added to the path of the script before running.   
          
      
