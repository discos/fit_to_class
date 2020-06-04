#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 12:15:50 2020

@author: debbio
"""
from fitslike_commons import Fitslike_commons


class Fitslike_component():
    """Fitslike semanti component interface"""

    def __init__(self, p_componentName, p_keywordsList):
        """
        Store fitslike keyword dictionary

        Parameters
        ----------
        p_ componentName : component name
            String id of this object
        p_keywordsList : keyword list
            List of this object fitslike keyword

        Returns
        -------
        None.

        """
        self.m_name = p_componentName
        self.m_helper = Fitslike_commons()
        self.m_keywordsDict = self.m_helper \
            .build_attribute_dict(p_keywordsList)

    def update(self, p_key, p_value):
        """
        Updates keyword value and type

        Parameters
        ----------
            p_key: keyword
            p_value: values

        Returns
        -------
        None.

        """
        if p_key in self.m_obsdata:
            self.m_keywordsDict[p_key] = p_value

    def dump(self):
        """Class usefull cotents dumping methods"""
        for l_key in self.m_keywordsDict.keys():
            self.m_helper.dump_attribute(l_key, self.m_keywordsDict[l_key])
