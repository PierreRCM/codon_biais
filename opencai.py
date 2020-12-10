import pandas as pd
import os
import matplotlib.pyplot as plt

df = pd.read_csv(os.getcwd()+"/CAI_EC", index_col=0).iloc[:, 0].astype(float).round(3)
# read CAI file