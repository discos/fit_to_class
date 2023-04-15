# fit_to_class

## Dependencies

Python deps:
 
 * astropy (>4)
 * numpy (==1.21)
 * traceback (opt)
 * memory_profiler 
    
To install deps:

      pip install -r requirements


##Installation 

    
      
      python setup.py install
      #or
      python setup.py develop  (to ve verified)
      
 ##Example usage 
 
      fit_to_class -s scan_conf/conf.json -p 12  -g 121ON,0OFF,0CAL  20210318-152823-3-20-CEPH_DEC/
      
      -s pipeline configuration file
      -p parallel processes working on file parsing
      -g geometry NumON,NumOFF,NumCAL, input files has to match the given order
      -fd feed to be parsed (optional, default all, to be verified)
      -r raw processing without grouping input file by geometry (legacy parameter, non really usefull now)
      -t scan type, e.g. OTF,ONOFF.. (legacy parmater not usefull now)
      
      Last argument is the relative folder path with input files
      
Test example usage:
      
    python test_fitslike_handler.py -i FitszillaFolder -p 4
    where:
        -i Input folder
        -p How many input files to be processed in parallel
  
    Output:
        - fit_to_class Folder inside Input Folder, usefull output data at InputFolder/fit_to_class/classfits/class_feed_x
    
Tested with: Python 3.8.5

 #Pipeline configuration file
 
The converter has an internal pipeline to solve the conversion, each stage might be disabled (mainly for development pourposes), the only case when the conf file comes in handy  is when we want to parse an OTF file where we want to avoid calibration when cal data aren't really availabe.
In this particular case just set to false the following entries:

    "data_calibration_counts":{
        "enabled": false
    },

    "data_calibration":{
        "enabled": false
    },
 
 Another usefull case might be the need of output spectrum measured in COUNTS:
 
    "classfit_conversion_counts":{
        "enabled": false
    },
 
 Summarizing, when we want to have output in COUNTS (with cal available):
 
    "data_calibration_counts":{
        "enabled": true
    },

    "classfit_conversion_counts":{
        "enabled": true
    },

on the other hand when we need the usual calibrated data:
   
    "data_calibration":{
        "enabled": true
    },

    "classfit_conversion":{
        "enabled": true
    }

## release notes

* 0.1 first commit
* 0.1.1 enhanced logging 


