import pdb
import json
from enum import Enum
from fit_to_class import fitslike_handler
from fit_to_class import fitslike_commons
import scanoptions


class PipelineTasks(Enum):
    """ Pipeline enum tasks """
    FITSLIKE_HANDLER= "fitslike_handler"
    SUBSCAN_PARSING= "subscan_parsing"
    GEOMETRY_GROUPING= "geomertry_grouping"
    SIGNAL_GROUPING= "signal_grouping"
    DATA_CALIBRATION= "data_calibration"
    CLASSFITS_CONVERSION ="classfit_conversion"

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
        _list.append(( PipelineTasks.FITSLIKE_HANDLER,
                    {     
                        "task": self._pipeline_fitslike_handler,
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
        # Geometry grouping
        _list.append((PipelineTasks.GEOMETRY_GROUPING,
                    {
                        "task": self._pipeline_geometry_grouping,
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
        _list.append((PipelineTasks.DATA_CALIBRATION,
                    {            
                        "task": self._pipeline_fitslike_handler,
                        "enabled" : False
                    })
        )
        # Classfits conversion
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
        # Fitsklike handler, mandatory
        _task_ft= self._pipeline_find_task(PipelineTasks.FITSLIKE_HANDLER)        
        if _task_ft:            
            _enabled= self._conf_task_is_enabled(PipelineTasks.FITSLIKE_HANDLER)
            if _enabled: 
                self._pipeline_task_set_enabled(_task_ft)            
            self._scan_pipeline.append(_task_ft)                        
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
        else:
            # Geometry grouping
            _task_geo= self._pipeline_find_task(PipelineTasks.GEOMETRY_GROUPING)
            if _task_geo:
                _enabled= self._conf_task_is_enabled(PipelineTasks.GEOMETRY_GROUPING)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_geo)
                self._scan_pipeline.append(_task_geo) 
            # Signal ref grouping
            _task_signal= self._pipeline_find_task(PipelineTasks.SIGNAL_GROUPING)
            if _task_signal:
                _enabled= self._conf_task_is_enabled(PipelineTasks.SIGNAL_GROUPING)
                if _enabled: 
                    self._pipeline_task_set_enabled(_task_signal)
                self._scan_pipeline.append(_task_signal) 
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
        self._logger.info("PIPELINE CLOSING\n")

    # PIPELINE TASKS

    def _pipeline_fitslike_handler(self, p_scan_context) -> None:
        """
        Creates fitlike handler 
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.FITSLIKE_HANDLER))
        _fh= fitslike_handler.Fitslike_handler( self._scan_options.raw,\
                    self._scan_options.geometry,\
                    self._scan_options.type,\
                    self._scan_options.feed,\
                    self._scan_options.parallel)   
        _fh.setOutputPath(self._scan_options.get_output_path())
        p_scan_context['fh']= _fh

    def _pipeline_scan_data(self, p_scan_context) -> None:
        """
        File scan operations
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.SUBSCAN_PARSING))
        try:
            _fh= p_scan_context['fh']
            _fh.scan_data(self._scan_options.folder)
        except Exception as e:
            self._logger.error(f"Exception on task execution {str(e)}")
    
    def _pipeline_geometry_grouping(self, p_scan_context) -> None:
        """
        Subscan group by geometry definitions
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.GEOMETRY_GROUPING))

    def _pipeline_signal_grouping(self, p_scan_context) -> None:
        """
        Subscan find on off and cal inside geometry groups or
        among the whole subscan list
        Also check if geometry definition fits subscans data, if 
        geometry is given.
    
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.SIGNAL_GROUPING))

    def _pipeline_normalize(self,p_scan_context) -> None:
        """
        Normalize data
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.DATA_CALIBRATION))
    
    def _pipeline_classfits(self,p_scan_context) -> None:
        """
        Classfits conversions
        Raise exception if errors

        p_scan_contex: dict, pipeline data
        """
        self._logger.info("PIPELINE EXECUTING: {}\n".format(PipelineTasks.CLASSFITS_CONVERSION))

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