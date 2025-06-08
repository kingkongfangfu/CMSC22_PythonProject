import numpy as np
import pandas as pd
import math

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import average, dendrogram
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import re
import os
import subprocess

def runGNKDST():
    global alleleFreqDF
    global dispanDF
    global dispanFile
    global namesOfLoci
    global numberOfLoci
    global numberOfIndividuals
    global numberOfPopulations
    global individualsPerPopn

    print("\nRunning DISPAN (GNKDST.exe)...")  # For demonstration, print a message to the console

    GNKDSTDirectory = os.path.join(os.path.dirname(__file__), 'DISPAN')
    executableGNKDST = os.path.join(GNKDSTDirectory, 'GNKDST')  # Adjust the name as needed

    if os.path.exists(executableGNKDST):
        try:
            inputFile = os.path.join(GNKDSTDirectory, "input")
            with open(inputFile, 'w', encoding='utf-8') as input_file:
                input_file.write(dispanFile)  # Write the opened structureFile content to input.txt
            
            # Change the working directory to ./structure/ and run the executable
            result = subprocess.run([executableGNKDST], cwd=GNKDSTDirectory, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Printing stdout and stderr for more info
            print("\nStructure executed successfully.")
            print("\nOutput:", result.stdout)
            print("\nError:", result.stderr)

        except subprocess.CalledProcessError as e:
            print(f"\nAn error occurred while running Structure: {e}")
            print("\nOutput:", e.stdout)
            print("\nError:", e.stderr)
            
        except Exception as e:
            print(f"\nAn error occurred: {e}")
    else:
        print("\nExecutable file not found.")
