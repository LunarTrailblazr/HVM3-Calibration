# David R Thompson
import numpy as np
import os, os.path, sys
import json


class HVM3Config:

  def __init__(self, filepath=None, mode='directsolar_wide'):

      if filepath is None:
          raise ValueError('Configuration file required')

      with open(filepath,'r') as fin:
          config_dict = json.load(fin)

          # Select the proper mode
          config_dict = config_dict[mode]
          for key, val in config_dict.items():
              setattr(self,key,val) 
      
      # Translate local filepaths to absolute where needed
      if hasattr(self,'base_directory'):
        my_directory = self.base_directory
      else:
        my_directory, my_executable = os.path.split(os.path.abspath(__file__))
      for fi in dir(self):
          if '_file' in fi:
              path = getattr(self,fi)
              if path[0] != '/':
                  path = os.path.join(my_directory, path)
                  setattr(self,fi,path)

      # Define masked columns
      if hasattr(self,'masked_columns'):
          ranges = [np.arange(a,b) for a,b in self.masked_columns] 
          self.masked_cols = np.concatenate(ranges, axis=0)
      else:
          self.masked_cols = []

