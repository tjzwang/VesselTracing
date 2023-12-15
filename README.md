# VesselTracing
  
Protocol:   

1. Obtain Image Coordinates for Data Analysis:

  File: FINAL_coordinate_exporter.m (or SHUFFFLED_FINAL_coordinate_exporter.m for shuffled data)
   
2. Calculate the distance of objects or pixels along the length of and away from blood vessels  
       Files:
   
       intensity_distance_analysis_FINAL.m (Use this file if data is correctly oriented. You can check orientation of data by looking at figures           that appear. In figure 4, the red vessel should be visible inside the blue tau data)

       object_distance_analysis_FINAL.m
       
    
3. Bin data along and away from vessel:
   
      Files:  
         For Intensity: intensity_binning_FINAL.m  
         For Objects: object_binning_FINAL.m
         To get volumes for object analysis: intensity_distance_analysis_for_volume.m

4. Get data for plots

     Files:
     INTENSITY_stat_plots_with_export.m
     OBJECT_stat_plots_with_export.m
  
NOTE: The folder "Necessary Functions" must be added to the path for scripts to run!
NOTE: All files must start in the format sampleID_Vessel_vesselID_. The file name can vary after the final underscore, ex. "Human2267_Vessel_3_DT". Failure to follow this format will result in errors, as sample and vessel IDs are extracted from the original file name. 
  
MATLAB Code Descriptions:

1.    FINAL_coordinate_exporter.m
      This script determines the X,Y,and Z position of each pixel in .tiff volumes and their corresponding values for inputed channels.                   Specifically in this case, it outputs the imaging data within 100 microns from the blood vessel surface of isolated blood vessels. This            must be run once per imaging sample. 

      Inputs:
          Data:
          1. Folder containing distance transforms of isolated objects. In other words, a folder containing images where pixel intensity                        corresponds to the distance from some isolated object. For example, distance transforms were created for isolated blood vessels,                   resulting in one distance transform image per vessel. These can be obtained in Imaris.

          2. .tif volume with pixel intensity corresponding to the cortical layer in which they reside. For examples, pixels in layer 1, should                 have an intensity of 1, then cortical layer 2 should have intensity of 2, etc. These images can be created in Imaris by manually                   drawing surfaces to encompass each cortical layer, then masking those surfaces and setting the pixel intensity within them to the                  layer number.

          3. .tif volume of the AT8 tau imaging channel, or the type of pathology you are trying to analyze relative to traced objects.

          ****Input folder names should be formatted such that the sample number is first, followed by an underscore. Coordinate files will be               saved as (SAMPLE #)_COORD

          Values:
          1. Save final data? Select how you would like the data to be outputted. Data can be saved as a .txt or .xlsx file, although .txt is                   recommended.

          2. Enter the folder in which you would like the coordinates to be saved. One file will be saved for each traced object. Ex. one                       coordinate file will be made for each isolated blood vessel. 
   
   Output Format:  
       Column: [  1                     2              3               4                  5                      6        ]  
               [  X(in pixels)    Y(in pixels)    Z(in pixels)     Intensity    DistanceFromSurface       Cortical Layer  ] 
  
2.    intensity_distance_analysis_FINAL.m    
      This script calculates the distance of pixels along and away from blood vessels. 
         
      Inputs:  
          Data:  
             1. Coordinates obtained in step 1.
             
          Values:  
             1. Image Dimensions in Microns (Includes lower and upper bounds for X, Y, and Z.)  
             2. Image Dimensions in Pixel Count (Include X dimension, Y dimension, and # of Z-stacks for each input: tau data, vessel data,                        and distance transform data)  
             3. Export path, name, and file type for exported data files. (.txt is recommended)  
               
      Output Format:  
       Column: [       1                 2               3            4                  5                      6                   7          ]  
               [  X(in microns)   Y(in microns)   Z(in microns)    Intensity     DistanceAlongVessel   DistanceFromSurface   Cortical Layer    ]   
       Data exports as .txt or .xlsx (.txt is recommended)  
               
           
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
          2. Folder containing files with calculated pixels distances of the distance transforms along and away from each blood vessel (Outputs from intensity_distance_analysis_for_volume.m)    
             
     Values:   
          1. Max Distance From Vessel (in microns)  
          2. Min Distance From Vessel (in microns)  
          3. Distance from Vessel Resolution (in microns)   
          4. Binning Along Vessel Frame Width (in microns) 
               
     Output Format:   
      Column: [    1        2            3                4                5-(Max Distance From Vessel - Min Distance From Vessel  ]  
              [ Vessel#  Frame#  MinFrameDistance    MaxFrameDistance           (Density For Each Distance From Vessel)            ]  
      Data exports as a .csv  


      
