import os
import logging
from fit_to_class import fitslike_handler
from fit_to_class import fitslike_commons
import scanoptions


class ScanPipeline:
    """
    Scan sequence pipeline
    """
    def __init__(self, p_logger) -> None:
        self._logger= p_logger
        self._scan_options= None

    def set_scan_options(self, p_options) -> None:
        """ 
        Setter scan options
        """
        if not isinstance(p_options, scanoptions.ScanOptions):
            self._logger.error("Wrong scan option instance type")
            self._scan_options= None
            return
        self._scan_options= p_options
        

    def scan(self) -> None:
        """
        Working pipeline         
        """        
        if self._scan_options.get_errors():
            self._logger.error("Errors on scan options parameters. Cannot start pipeline.")
            return
        # Output path building
        l_fh= fitslike_handler.Fitslike_handler( self._scan_options.raw,\
                    self._scan_options.type,\
                    self._scan_options.feed,\
                    self._scan_options.parallel)        
        
        l_fh.setOutputPath(self._scan_options.get_output_path())
        # Data scan
        l_fh.scan_data(self._scan_options.folder)
        if self._scan_options.raw:
        # Data conversion
            l_fh.group_on_off_cal()
            l_fh.ClassFitsAdaptations()
        else:
            l_fh.group_on_off_cal()
            l_fh.normalize()
            l_fh.ClassFitsAdaptations()