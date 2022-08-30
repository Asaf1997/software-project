import numpy as np
import pandas as pd

a = pd.read_csv('spk_0.txt', header=None)
print(len(a))
print()
print(len(a[0]))


c = pd.read_csv('test_transcript_c.txt', header=None)
print(c)