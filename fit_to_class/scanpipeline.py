import pdb
import os
import sys
import shutil
import traceback
import inspect
import json
import numpy as np
from enum import Enum
from fit_to_class import fitslike_handler
from fit_to_class import fitslike_commons
from fit_to_class.fitslike_scantype import ScanType
from fit_to_class.fitslike_scantype import ScanTypes
from fit_to_class import scanoptions

class PipelineTasks(Enum):
    """ Pipeline enum tasks """
    SUBSCAN_LIST= "subscan_list"
    PRE_PARSING= "pre_parsing"
    FITSLIKE_HANDLER= "fitslike_handler"
    SUBSCAN_PARSING= "subscan_parsing"
    GEOMETRY_GROUPING= "geomertry_grouping"
    SIGNAL_GROUPING= "signal_grouping"
    DATA_CALIBRATION= "data_calibration"
    DATA_CALIBRATION_COUNTS= "data_calibration_counts"
    CLASSFITS_CONVERSION ="classfit_conversion"
    CLASSFITS_CONVERSION_COUNTS ="classfit_conversion_counts"

class ScanPipeline:
    """
    Scan sequence pipeline
    """
    # ENUM
    FITSLIKE_HANDLER= "fitslike_handler"

    # PUBLIC

    def __init__(self, p_logger) -> None:
        # Logger
        self._logger= p_logger
        # Input parameters
        self._scan_options= None
        self._conf_json= {}
        # Pipeline conf and data
        self._pipeline_task_list= self.pipeline_init()        
        self._scan_pipeline= []
        self._scan_context= {}
        self._pipeline_errors= False
        self._scan_list= {
            'summary': '',
            'subscan': []
        }
        self._scan_list_group= []


    def set_scan_options(self, p_options) -> None:
        """ 
        Setter scan options
        """
        if not isinstance(p_options, scanoptions.ScanOptions):
            self._logger.error("Wrong scan option instance type")
            self._scan_options= None
            return
        self._scan_options= p_options
        # Read and parse scan conf file
        self._conf_read()        
        # Build pipeline
        self._pipeline_build()
        # Print pipeline 
        self._pipeline_print()

    def pipeline_start(self) -> None:
        """
        Working pipeline         
        """        
        if self._scan_options.get_errors():
            self._logger.error("Errors on scan options parameters. Cannot start pipeline.")
            return
        # TODO review pipeline tasks with new methods
        self._pipeline_errors= False
        for task in self._scan_pipeline:
            # check for critical errors
            if self._pipeline_errors:
                self._logger.info("Stopping pipe line execution due to critical errors")
                self._pipeline_close()
                return
            # Task reading
            task_enabled= self._pipeline_task_is_enabled(task)
            #self._logger.info("enabled:{} - {}".format(self._task_is_enabled(task), self._pipeline_task_name(task)))
            if not task_enabled:
                continue
            # Task execution
            self._pipeline_task_execute(task)
        # Closing
        #self._logger.info("Closing pipeline")
        self._pipeline_close()             

    # BUILD PIPELINE

    def pipeline_init(self) -> list:
        """
        Pipeline task inital set up
        [(name, task dict),(name, task dict),..]
        """
        _list=[]
        self._scan_list= {
            'summary': '',
            'subscan': []
        }
        self._scan_list_group= []
        # Filter subscan folder
        _list.append((PipelineTasks.SUBSCAN_LIST,
                    {     
                        "task": self._pipeline_subscan_list,
                        "enabled" : False
                    })
        )
        # Pre parsing summary
        _list.append((PipelineTasks.PRE_PARSING,
                    {     
                        "task": self._pipeline_pre_parsing,
                        "enabled" : False
                    })
        )
        # Geometry grouping
        _list.append((PipelineTasks.GEOMETRY_GROUPING,
                    {
                        "task": self._pipeline_geometry_grouping,
                        "enabled" : False
                    })
        )
        # File scan
        _list.append(( PipelineTasks.SUBSCAN_PARSING,
                    {                    
                        "task": self._pipeline_scan_data,
                        "enabled" : False
                    })
        )        
        # Signal reference grouping        
        _list.append((PipelineTasks.SIGNAL_GROUPING,
                    {         
                        "task": self._pipeline_signal_grouping,
                        "enabled" : False
                    })
        )
        # Normalize values 
        _list.append((PipelineTasks.DATA_CALIBRATION_COUNTS,
                    {            
                        "task": self._pipeline_normalize_counts,
                        "enabled" : False
                    })
        )
        # Classfits conversion counts
        _list.append((PipelineTasks.CLASSFITS_CONVERSION_COUNTS,
                    {                    
                        "task": self._pipeline_classfits_counts,
                        "enabled" : False
                    })
        )
        # Normalize values
        _list.append((PipelineTasks.DATA_CALIBRATION,
                    {            
                        "task": self._pipeline_normalize,
                        "enabled" : False
                    })
        )

        # Classfits conversion full
        _list.append((PipelineTasks.CLASSFITS_CONVERSION,
                    {                    
                        "task": self._pipeline_classfits,
                        "enabled" : False
                    })
        )        
        return _list

    def _pipeline_find_task(self, p_name) ->dict:
        """
        Task look up from all task list
        """
        _task_name= p_name
        #pdb.set_trace()
        for tup in self._pipeline_task_list:
            if _task_name == tup[0]:
                return tup[1]
        return None

    def _pipeline_build(self) -> None:
        """
        Pipeline building from scan options
        and from scan conf file
        
        # # Output path building
        # l_fh= fitslike_handler.Fitslike_handler( self._scan_options.raw,\
        #             self._scan_options.geometry,\
        #             self._scan_options.type,\
        #             self._scan_options.feed,\
        #             self._scan_options.parallel)        
        
        # l_fh.setOutputPath(self._scan_options.get_output_path())
        # # Data scan
        # l_fh.scan_data(self._scan_options.folder)
        # # @TODO Check scan pipeline in for OTF/maps
        # if self._scan_options.raw:
        # # Data conversion
        #     l_fh.group_on_off_cal()
        #     l_fh.ClassFitsAdaptations()
        # else:
        #     l_fh.group_on_off_cal()
        #     l_fh.normalize()
        #     l_fh.ClassFitsAdaptations()

        """
        self._scan_pipeline=[]        
        # Subscan folder filter, mandatory
        _task_sublist= self._pipeline_find_task(PipelineTasks.SUBSCAN_LIST)        
        if _task_sublist:            
            _enabled= self._conf_task_is_enabled(PipelineTasks.SUBSCAN_LIST)
            if _enabled: 
                self._pipeline_task_set_enabled(_task_sublist)            
            self._scan_pipeline.append(_task_sublist)    
        # Summary pre scan
        _task_pre_parsing= self._pipeline_find_task(PipelineTasks.PRE_PARSING) 
        if _task_pre_parsing:            
            _enabled= self._conf_task_is_enabled(PipelineTasks.PRE_PARSING)
            if _enabled: 
                self._pipeline_task_set_enabled(_task_pre_parsing)
            self._scan_pipeline.append(_task_pre_parsing)    
        # Geometry grouping        
        _task_geo= self._pipeline_find_task(PipelineTasks.GEOMETRY_GROUPING)
        if _task_geo:
            _enabled= self._conf_task_is_enabled(PipelineTasks.GEOMETRY_GROUPING)
            if _enabled: 
                self._pipeline_task_set_enabled(_task_geo)
            self._scan_pipeline.append(_task_geo)                        
        # Scan files, mandatory
        _task_scan= self._pipeline_find_task(PipelineTasks.SUBSCAN_PARSING)        
        if _task_scan:
            _enabled= self._conf_task_is_enabled(PipelineTasks.SUBSCAN_PARSING)
            if _enabled: 
                self._pipeline_task_set_enabled(_task_scan)
            self._scan_pipeline.append(_task_scan) 
        # First branch, raw
        if self._scan_options.raw:
            # Signal ref grouping
            _task_signal= self._pipeline_find_task(PipelineTasks.SIGNAL_GROUPING)
            if _task_signal:
                _enabled= self._conf_task_is_enabled(PipelineTasks.SIGNAL_GROUPING)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_signal)
                self._scan_pipeline.append(_task_signal) 
            # Classfits conversion
            _task_classfits= self._pipeline_find_task(PipelineTasks.CLASSFITS_CONVERSION)
            if _task_classfits:
                _enabled= self._conf_task_is_enabled(PipelineTasks.CLASSFITS_CONVERSION)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_classfits)
                self._scan_pipeline.append(_task_classfits) 
            # Classfits conversion counts
            _task_classfits_counts= self._pipeline_find_task(PipelineTasks.CLASSFITS_CONVERSION_COUNTS)
            if _task_classfits_counts:
                _enabled= self._conf_task_is_enabled(PipelineTasks.CLASSFITS_CONVERSION_COUNTS)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_classfits_counts)
                self._scan_pipeline.append(_task_classfits_counts) 
        else:             
            # Signal ref grouping
            _task_signal= self._pipeline_find_task(PipelineTasks.SIGNAL_GROUPING)
            if _task_signal:
                _enabled= self._conf_task_is_enabled(PipelineTasks.SIGNAL_GROUPING)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_signal)
                self._scan_pipeline.append(_task_signal) 
            # Data normalization counts
            _task_norm_counts= self._pipeline_find_task(PipelineTasks.DATA_CALIBRATION_COUNTS)
            if _task_norm_counts:
                _enabled= self._conf_task_is_enabled(PipelineTasks.DATA_CALIBRATION_COUNTS)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_norm_counts)
                self._scan_pipeline.append(_task_norm_counts) 
            # Classfits conversion counts
            _task_classfits_counts= self._pipeline_find_task(PipelineTasks.CLASSFITS_CONVERSION_COUNTS)
            if _task_classfits_counts:
                _enabled= self._conf_task_is_enabled(PipelineTasks.CLASSFITS_CONVERSION_COUNTS)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_classfits_counts)
                self._scan_pipeline.append(_task_classfits_counts) 
            # Data normalization
            _task_norm= self._pipeline_find_task(PipelineTasks.DATA_CALIBRATION)
            if _task_norm:
                _enabled= self._conf_task_is_enabled(PipelineTasks.DATA_CALIBRATION)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_norm)
                self._scan_pipeline.append(_task_norm) 
            # Classfits conversion
            _task_classfits= self._pipeline_find_task(PipelineTasks.CLASSFITS_CONVERSION)
            if _task_classfits:
                _enabled= self._conf_task_is_enabled(PipelineTasks.CLASSFITS_CONVERSION)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_classfits)
                self._scan_pipeline.append(_task_classfits) 
            
        
    def _pipeline_task_is_enabled(self, p_task) -> bool:
        """
        Check if given task is enabled

        p_task: dict, pipeline task
        """
        _key= "enabled"
        if type(p_task) is not dict:
            return False
        if _key in p_task:
            return p_task[_key]
        return False

    def _pipeline_task_set_enabled(self, p_task) -> None:
        """
        Setter pipeline task enabled
        """
        _key= "enabled"
        if type(p_task) is not dict:
            raise TypeError("Given task not a dict")            
        if _key in p_task:            
            p_task[_key]= True
        else:
            raise KeyError("{} not present on task {}".format(_key,json.dumps(p_task)))            

    def _pipeline_task_name(self, p_task) -> str:
        """
        Returns task name
        """    
        _key= "task"
        if type(p_task) is not dict:
            return False
        if _key in p_task:
            return p_task[_key].__name__
        return "unknown"

    def _pipeline_task_execute(self, p_task) -> None:
        _key= "task"
        if type(p_task) is not dict:
            self._logger.error("Given task in not a dict!")
            self._pipeline_errors= True
            return False
        if _key in p_task:
            try:
                p_task[_key](self._scan_context)
            except Exception as e:                                
                self._pipeline_errors= True
                self._logger.error("Error on task {}:\n {}".format(self._pipeline_task_name(p_task), str(e)))
                _exc_info= sys.exc_info()
                traceback.print_exc(*_exc_info)                
            finally:                
                return
        # Unknwon task
        self._logger.error("Cannot execute task {}".format(self._pipeline_task_name(p_task)))
        self._pipeline_errors= True

    def _pipeline_close(self):
        """
        Cleaning pipeline operations 
        """
        self._logger.info("PIPELINE CLOSING\n")

    # PIPELINE TASKS

    def _pipeline_subscan_list(self, p_scan_context) -> None:
        """
        Scan list first file look up
        It prepares master scan list
        _scan_list= {
            'summary': summary_file
            'subscan': [scubscan file list]
        }
        """
        folder= self._scan_options.folder
        files_fits= [os.path.join(folder,f) for f in os.listdir(folder) if f.endswith('.fits')]
        self._scan_list['summary']= [f for f in files_fits if '/sum' in f.lower()]
        # Summary check if exists (it has to be one only)
        if  not self._scan_list['summary']:
            self._logger.error(f"Missing summary.fits from {folder}")
            self._pipeline_errors= True
            return
        # Summary, from list to single file
        self._scan_list['summary']= self._scan_list['summary'][0]
        # Subscans
        self._scan_list['subscan']= [f for f in files_fits if 'summary' not in f]
        # Print
        #self._subscan_list_print()
            
    def _pipeline_pre_parsing(self, p_scan_context) -> None:
        """
        Before joining subscans by geo it looks into summary.fits
        for:
            - geometry keyword
            - scan type keyword
        """
        self._logger.info("\nPIPELINE EXECUTING: {}".format(PipelineTasks.PRE_PARSING))
        # Check summary enlisted
        if 'summary' not in self._scan_list:
            self._logger.error(f"Missing summary.fits from input folder!")            
            raise Exception("Missing summary.fits from input folder!")
        # scanning only summary
        if self._scan_list['summary']:
            self._logger.info(f"Pre parsing summary.fits")                                    
            _fh= fitslike_handler.Fitslike_handler( self._scan_options.raw,\
                                                    self._scan_options.geometry,\
                                                    self._scan_options.type,\
                                                    self._scan_options.feed,\
                                                    self._scan_options.parallel)                        
            _summary_group= {}
            _summary_group['summary']= self._scan_list['summary']
            _summary_group['subscan']= []
            try:
                _fh.scan_data(_summary_group)
            except Exception as e:                
                self._logger.error(f"Exception on summary pre scan {str(e)}")
                raise Exception(f"Exception on summary pre scan {str(e)}")
            # Review scan options (geometry , type)
            # Geometry from command line has priority over geometry from summary                        
            _geo_summary= _fh.get_geometry_definition()            
            self._logger.info(f"Geometry definition from summary.fits: {_geo_summary}")
            if not self._scan_options.geometry.get_geometry():
                # Check if summary has geometry info
                if not _geo_summary:
                    raise Exception("Missing geometry input definition!")
                self._scan_options.geometry= _geo_summary
                self._logger.info("Using geometry definition from summary.fits")                
            else: 
                self._logger.info("Using geometry definition from command line")                                    
            # Scan type
            _scan_type_summary= _fh.get_scan_type()
            self._logger.info(f"Scan type definition from summary.fits: {_scan_type_summary}")
            if not ScanType.is_valid(self._scan_options.type):
                # Check if summary has scan_type info
                if not _scan_type_summary:
                    raise Exception("Missing scan type input definition!")
                self._scan_options.type= _scan_type_summary
                self._logger.info("Using scan type definition from summary.fits")
            else: 
                self._logger.info("Using scan type definition from command line")                        
            

    def _pipeline_geometry_grouping(self, p_scan_context) -> None:
        """        
        Subscan group by geometry definitions
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("\nPIPELINE EXECUTING: {}".format(PipelineTasks.GEOMETRY_GROUPING))
        # TODO Leggere i file dalla folder e sapararli secondo geometria,
        # Dovrebbe essere possibile passare i singoli gruppi allo scan di FH 
        _subcan_len= len(self._scan_list['subscan'])
        _geo_loop_len= self._scan_options.geometry.get_loop_len()
        _mod= _subcan_len % _geo_loop_len        
        if _mod != 0:
            self._logger.error(f"Unexpected susbcan number, cannot gorup by given geometry:\n"\
                + f"Total subscan number: {_subcan_len}\n"\
                + f"Geometry single group len: {_geo_loop_len}\n"\
                + f"MOD is not zero: {_mod}")
            self._pipeline_errors= True
            return
        # Ok subscan files number fits geo length        
        _geo_group_num= int(_subcan_len / _geo_loop_len)
        self._scan_list['subscan']= [g.tolist() for g in np.array_split(self._scan_list['subscan'], _geo_group_num)]
        # Creating self contained scan group list
        _cnt=0
        for el in self._scan_list['subscan']:
            if type(el) is list:
                self._scan_list_group.append({
                    'group_id': _cnt,
                    'summary': self._scan_list['summary'],
                    'subscan': el
                })
                _cnt= _cnt +1            
        if _cnt == 0:
            # No geo, one single group
            self._scan_list_group.append({
                    'group_id': _cnt,
                    'summary': self._scan_list['summary'],
                    'subscan': self._scan_list['subscan']
                })
        # Print geo groups
        #self._subscan_list_print()
        self._subscan_group_print()
   
    
    def _pipeline_scan_data(self, p_scan_context) -> None:
        """
        File scan operations
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """        
        self._logger.info("\nPIPELINE EXECUTING: {}".format(PipelineTasks.SUBSCAN_PARSING))
        # Splitting working load based on geo sub list
        # creating a fitslike_handler scan for every geo group
        if self._scan_list_group:
            # Geometry applied
            for _geo_group in self._scan_list_group:                
                self._logger.info(f"Ready to scan group {_geo_group['group_id']}")
                # Handler for every geo group
                _geo_group['has_errors'] = False
                _fh= fitslike_handler.Fitslike_handler( self._scan_options.raw,\
                                                        self._scan_options.geometry,\
                                                        self._scan_options.type,\
                                                        self._scan_options.feed,\
                                                        self._scan_options.parallel)                                                                        
                _fh.setOutputPath(self._scan_options.get_output_path().rstrip('/')+'_'+str(_geo_group['group_id']))
                _geo_group['fh']= _fh
                try:
                    _fh.scan_data(_geo_group)
                except Exception as e:                    
                    self._logger.error(f"Exception on task execution {str(e)}")
                    #_exc_info= sys.exc_info()                    
                    #traceback.print_exc()
                    _geo_group['has_errors']= True       

        
    def _pipeline_signal_grouping(self, p_scan_context) -> None:
        """
        Subscan find on off and cal inside geometry groups or
        among the whole subscan list
        Also check if geometry definition fits subscans data, if 
        geometry is given.
    
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("\nPIPELINE EXECUTING: {}".format(PipelineTasks.SIGNAL_GROUPING))
        # Traversing scan groups
        if self._scan_list_group:            
            for _geo_group in self._scan_list_group:             
                self._logger.info(f"On Off Call for group {_geo_group['group_id']}")                   
                # FH, fitslike_handler contains scan data for the given geo group
                if 'fh' in _geo_group:
                    try:
                        _geo_group['fh'].group_on_off_cal()
                    except Exception as e:
                        self._logger.error(f"Exception on data during On Off Cal step {str(e)}")
                        self._logger.error("\nEXCEPTION\n-------------------------------------")
                        _exc_info= sys.exc_info()                    
                        traceback.print_exc()
                        self._logger.error("\n-------------------------------------")
                else:
                    self._logger.error(f"Missing fitslike handler for group {_geo_group['group_id']}")                    


    def _pipeline_normalize_counts(self,p_scan_context) -> None:
        """
        Normalize data
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.DATA_CALIBRATION_COUNTS))
        # Traversing scan groups
        if self._scan_list_group:            
            for _geo_group in self._scan_list_group:             
                self._logger.info(f"On Off Call for group {_geo_group['group_id']}")                   
                # FH, fitslike_handler contains scan data for the given geo group
                if 'fh' in _geo_group:
                    try:
                        _geo_group['fh'].normalize("counts")
                    except Exception as e:
                        self._logger.error(f"Exception on data counts {str(e)}")
                        self._logger.error("\nEXCEPTION\n-------------------------------------")
                        _exc_info= sys.exc_info()                    
                        traceback.print_exc()
                        self._logger.error("\n-------------------------------------")
                else:
                    self._logger.error(f"Missing fitslike handler for group {_geo_group['group_id']}")            

    def _pipeline_normalize(self,p_scan_context) -> None:
        """
        Normalize data
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.DATA_CALIBRATION))
        # Traversing scan groups
        if self._scan_list_group:            
            for _geo_group in self._scan_list_group:             
                self._logger.info(f"On Off Call for group {_geo_group['group_id']}")                   
                # FH, fitslike_handler contains scan data for the given geo group
                if 'fh' in _geo_group:
                    try:
                        _geo_group['fh'].normalize("cal")
                    except Exception as e:
                        self._logger.error(f"Exception on data calibration {str(e)}")
                        self._logger.error("\nEXCEPTION\n-------------------------------------")
                        _exc_info= sys.exc_info()                    
                        traceback.print_exc()
                        self._logger.error("\n-------------------------------------")
                else:
                    self._logger.error(f"Missing fitslike handler for group {_geo_group['group_id']}")    
        

    def _pipeline_classfits_counts(self,p_scan_context) -> None:
        """
        Classfits conversions from counts instead of Kelvin
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.CLASSFITS_CONVERSION_COUNTS))
        # Traversing scan groups
        if self._scan_list_group:            
            for _geo_group in self._scan_list_group:             
                self._logger.info(f"On Off Call for group {_geo_group['group_id']}")                   
                # FH, fitslike_handler contains scan data for the given geo group                
                if 'fh' in _geo_group:
                    try:
                        _geo_group['fh'].ClassFitsAdaptations("counts")
                    except Exception as e:
                        self._logger.error(f"Exception on data conversion to classfit counts {str(e)}")
                        self._logger.error("\nEXCEPTION\n-------------------------------------")
                        _exc_info= sys.exc_info()                    
                        traceback.print_exc()
                        self._logger.error("\n-------------------------------------")
                else:
                    self._logger.error(f"Missing fitslike handler for group {_geo_group['group_id']}")
    
    def _pipeline_classfits(self,p_scan_context) -> None:
        """
        Classfits conversions
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.CLASSFITS_CONVERSION))
        # Traversing scan groups
        if self._scan_list_group:            
            for _geo_group in self._scan_list_group:             
                self._logger.info(f"On Off Call for group {_geo_group['group_id']}")                   
                # FH, fitslike_handler contains scan data for the given geo group
                if 'fh' in _geo_group:
                    try:
                        _geo_group['fh'].ClassFitsAdaptations("cal")
                    except Exception as e:
                        self._logger.error(f"Exception on data conversion to classfit {str(e)}")
                        self._logger.error("\nEXCEPTION\n-------------------------------------")
                        _exc_info= sys.exc_info()                    
                        traceback.print_exc()
                        self._logger.error("\n-------------------------------------")
                else:
                    self._logger.error(f"Missing fitslike handler for group {_geo_group['group_id']}")

    # UTY 
                        
    def _conf_read(self):
        """
        Reading json pipeline conf file
        {
            "taskname":{
                "enabled" : true / false
            },
            ...
        }
        """        
        try:
            conf_file= open(self._scan_options.scan_conf_path,)
            self._conf_json= json.loads(conf_file.read())
            conf_file.close()
        except Exception as e:
            self._logger.error("Exception reading pipeline conf file: {} \n {}".\
                format(self._scan_options.scan_conf_path, str(e)))
            self._pipeline_errors= True
        finally:
            self._logger.info(f"\nPipeline task configuration:\n {json.dumps(self._conf_json, indent= 2, separators=(',', ':'))}\n")
            
    def _conf_parse(self):
        """
        Parsing scan conf dictionary
        Enable / disable pipeline steps
        """
        pass

    def _conf_task_is_enabled(self, p_task_name) -> bool:
        """
        Check if given taskname is enabled
        """
        _task_name= p_task_name.value
        _key= "enabled"        
        if _task_name in self._conf_json:
            if _key in self._conf_json[_task_name]:
                return self._conf_json[_task_name][_key]
        return False
    
    def _task_is_enabled(self, p_task) -> bool:
        """
        Return if the pipleine task is enabled
        """
        if type(p_task) is not dict:
            self._pipeline_errors= True
            self._logger.error(f"Given object is not a pipeline task dict!")
            return  False
        _key= "enabled"
        if _key not in p_task:
            self._pipeline_errors= True
            self._logger.error(f"Given pipeline task enabled key missing!")
            return False
        return p_task[_key]

    def _pipeline_print(self) -> None:
        """ Print pipeline """
        _active_task_cnt= 0        
        _body=""        
        for task in self._scan_pipeline:
            if self._task_is_enabled(task):
                _body = _body + f"{task['task'].__name__}\n"
                _active_task_cnt= _active_task_cnt +1
        _title= f"\nPIPELINE EXECUTION LIST ({_active_task_cnt} items):\n"
        self._logger.info(_title + _body)

    def _subscan_list_print(self) -> None:
        """ Print input file list """
        # Print file list        
        self._logger.info("SUSBCAN FILE LIST:")
        self._logger.info(f"Summary: {self._scan_list['summary']}")
        _cnt= 0 
        for f in self._scan_list['subscan']:
            if type(f) is list:
                self._logger.info("Geometry group:")
                for ff in f:
                    _cnt= _cnt +1            
                    self._logger.debug(f" Subscan {_cnt}: {ff}")
            else:
                _cnt= _cnt +1            
                self._logger.debug(f"Subscan {_cnt}: {f}")            
        self._logger.info(f"Subscan files : {_cnt}")

    def _subscan_group_print(self) -> None:
        """ Print scan group working list """
        # Print file list        
        self._logger.info("SUSBCAN GROUPS:")                
        _cnt= 0
        for _group in self._scan_list_group:
            self._logger.info(f"Group ID: {_group['group_id']}")
            _head, _tail= os.path.split(_group['summary'])
            self._logger.info(f"Summary: {_tail}")
            for _subscan in _group['subscan']:
                _head, _tail= os.path.split(_subscan)
                self._logger.debug(f"Subscan {_cnt}: {_tail}")            
                _cnt= _cnt +1            
        self._logger.info(f"Subscan group files number : {_cnt}")