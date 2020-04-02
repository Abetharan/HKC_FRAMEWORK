import os 
import string
import yaml
class Input:

    def __init__(self, path):
        with open(path, 'r') as stream:
            self.yaml_file = yaml.safe_load(stream)
