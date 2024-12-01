import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import seaborn as sns
import glob
import cv2
import uproot
import os

class pos_scan:
    def __init__(self, root_file_path):
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        file = uproot.open(os.path.join(self.script_dir, root_file_path))
        self.tree = file["tree"].arrays(library="np")

    def convert_x(self, x, x_min, width = 16):
        return (x - x_min)/width + 1/2

    def convert_y(self, y, y_max, width = 18):
        return (y_max - y)/width + 1/2

class pair:
    def __init__(self, val, err):
        self.val = val
        self.err = err