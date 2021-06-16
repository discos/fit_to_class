import sys
import json
from typing import List
from fit_to_class import fitslike_handler
from fit_to_class import fitslike_commons
import scanoptions


class ScanPipeline:
    """
    Scan sequence pipeline
    """

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
        self._read_scan_conf()
        self._parse_scan_conf()
        # Build pipeline
        self._pipeline_build()
        # Print pipeline 
        self._pipeline_print()

    def scan(self) -> None:
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
            self._logger.info("Pipeline::{} enabled:{}".format(self._pipeline_task_name(task)))
            if not task_enabled:
                continue
            # Task execution
            self._pipeline_task_execute(task)
        # Closing
        self._logger.info("Closing pipeline")
        self._pipeline_close()             

    # BUILD PIPELINE

    def pipeline_init(self) -> list:
        """
        Pipeline task inital set up
        [(name, task dict),(name, task dict),..]
        """
        _list=[]
        # Fitslike handler 
        _list.append(( "fitslike",
                    {     
                        "task": self._pipeline_fitslike_handler,
                        "enabled" : False
                    })
        )
        # File scan
        _list.append(( "scan",
                    {                    
                        "task": self._pipeline_scan_data,
                        "enabled" : False
                    })
        )
        # Geometry grouping
        _list.append(("geometry",
                    {
                        "task": self._pipeline_geometry_grouping,
                        "enabled" : False
                    })
        )
        # Signal reference grouping        
        _list.append(("on_off",
                    {         
                        "task": self._pipeline_signal_grouping,
                        "enabled" : False
                    })
        )
        # Normalize values
        _list.append(("normalize",
                    {            
                        "task": self._pipeline_fitslike_handler,
                        "enabled" : False
                    })
        )
        # Classfits conversion
        _list.append(("classfits",
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
        _task_name= str(p_name)
        for tup in self._pipeline_task_list:
            if _task_name in tup[0]:
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
        # Fitsklike handler, mandatory
        _task_ft= self._pipeline_find_task("fitliske")
        if _task_ft:            
            _enabled= self._conf_task_is_enabled("fitslike")
            if _enabled: 
                self._pipeline_task_set_enabled(_task_ft)
            self._scan_pipeline.append(_task_ft)            
        # Scan files, mandatory
        _task_scan= self._pipeline_find_task("scan")        
        if _task_scan:
            _enabled= self._conf_task_is_enabled("scan")
            if _enabled: 
                self._pipeline_task_set_enabled(_task_scan)
            self._scan_pipeline.append(_task_scan) 
        # First branch, raw
        if self._scan_options.raw:
            # Signal ref grouping
            _task_signal= self._pipeline_find_task("on_off")
            if _task_signal:
                _enabled= self._conf_task_is_enabled("on_off")
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_signal)
                self._scan_pipeline.append(_task_signal) 
            # Classfits conversion
            _task_classfits= self._pipeline_find_task("classifits")
            if _task_classfits:
                _enabled= self._conf_task_is_enabled("classfits")
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_classfits)
                self._scan_pipeline.append(_task_classfits) 
        else:
            # Geometry grouping
            _task_geo= self._pipeline_find_task("geometry")
            if _task_geo:
                _enabled= self._conf_task_is_enabled("geometry")
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_geo)
                self._scan_pipeline.append(_task_geo) 
            # Signal ref grouping
            _task_signal= self._pipeline_find_task("on_off")
            if _task_signal:
                _enabled= self._conf_task_is_enabled("on_off")
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_signal)
                self._scan_pipeline.append(_task_signal) 
            # Data normalization
            _task_norm= self._pipeline_find_task("normalize")
            if _task_norm:
                _enabled= self._conf_task_is_enabled("normalize")
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_norm)
                self._scan_pipeline.append(_task_norm) 
            # Classfits conversion
            _task_classfits= self._pipeline_find_task("classifits")
            if _task_classfits:
                _enabled= self._conf_task_is_enabled("classfits")
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
            return False
        if _key in p_task:
            try:
                p_task[_key]()
            except Exception as e:                
                self._logger.error("Error on task {}:\n {}".format(self._pipeline_task_name(p_task), str(e)))
                self._pipeline_errors= True
            finally:
                return
        # Unknwon task
        self._logger.error("Cannot execute task {}".format(self._pipeline_task_name(p_task)))
        self._pipeline_errors= True

    def _pipeline_close(self):
        """
        Cleaning pipeline operations 
        """

    # PIPELINE TASKS

    def _pipeline_fitslike_handler(self, p_scan_context) -> None:
        """
        Creates fitlike handler 
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """


    def _pipeline_scan_data(self, p_scan_context) -> None:
        """
        File scan operations
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """

    
    def _pipeline_geometry_grouping(self, p_scan_context) -> None:
        """
        Subscan group by geometry definitions
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """

    def _pipeline_signal_grouping(self, p_scan_context) -> None:
        """
        Subscan find on off and cal inside geometry groups or
        among the whole subscan list
        Also check if geometry definition fits subscans data, if 
        geometry is given.
    
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """

    def _pipeline_normalize(self,p_scan_context) -> None:
        """
        Normalize data
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
    
    def _pipeline_classfits(self,p_scan_context) -> None:
        """
        Classfits conversions
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """

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
            self._conf_json= json.loads(conf_file)
            conf_file.close()
        except Exception as e:
            self._logger.error("Exceptin reading pipeline conf file:\n {}".format(str(e)))

    def _conf_parse(self):
        """
        Parsing scan conf dictionary
        Enable / disable pipeline steps
        """

    def _conf_task_is_enabled(self, p_task_name) -> bool:
        """
        Check if given taskname is enabled
        """
        _task_name= str(p_task_name)
        _key= "enabled"
        if _task_name in self._conf_json:
            if _key in self.conf_json[_task_name]:
                return self.conf_json[_task_name][_key]
        return False

    def _pipeline_print(self) -> None:
        """ Print pipeline """
        print(json.dumps(self._scan_pipeline))