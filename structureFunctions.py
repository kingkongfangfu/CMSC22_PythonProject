import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import os
import subprocess

def openStructureFile():
    global structureFile
    global structureDataset
    global namesOfLoci
    global numberOfLoci
    global numberOfIndividuals
    global numberOfPopulations
    global individualsPerPopn
    
    file_path = filedialog.askopenfilename(initialdir=".", filetypes=[("Text files", "*.txt"), ("All files", "*.*")], title="Select a file")
    if file_path:
        with open(file_path, 'r', encoding='utf-8') as file:  # Ensure UTF-8 encoding for reading
            structureFile = file.read()
            if structureFile:  # Check if the file has content
                structureDataset = structureFile.splitlines()[1:]
                print(structureFile)  # For demonstration, print the content to the console
                namesOfLoci = structureFile.splitlines()[0].split()  # Extract the first line for locus names
                print("\nNames of loci:", namesOfLoci)  # Print the names of loci
                numberOfLoci = len(structureFile.splitlines()[0].split())
                print("\nNumber of loci:", numberOfLoci) 
                numberOfIndividuals = len(structureFile.splitlines()) // 2  # Assuming two lines per individual
                print("\nNumber of individuals:", numberOfIndividuals)
                uniquePopulations = sorted(set(line.split()[1] for line in structureDataset))  # Extract unique populations
                print("\nUnique populations:", uniquePopulations)  # Print unique populations
                numberOfPopulations = len(uniquePopulations)
                print("\nNumber of populations:", numberOfPopulations)  # Print the number of populations

                from collections import defaultdict
                expectedAlleles = numberOfLoci # expected number of allele 
                expectedStrings = 2 + expectedAlleles # expected number of strings per individual
                individualsPerPopn = defaultdict(set)  # population -> set of individuals
                for line in structureDataset:
                    if len(line.split()) != expectedStrings:
                        print(f"\nWarning: Line does not match expected format: {line.strip()}")
                        break # Stop processing if the line does not match the expected format
                    parts = line.split()
                    individual = parts[0]  # First part is the individual name
                    population = parts[1]  # Second part is the population name
                    individualsPerPopn[population].add(individual)  # Add individual to the set for the population
                
                # Print the counts
                for population, individuals in individualsPerPopn.items():
                    print(f"\nPopulation {population} has {len(individuals)} unique individuals.")

            else:
                print("\nThe file is empty.")  # Inform if the file is empty

def runStructure():
    # Placeholder for the runStructure function
    print("\nRunning Structure...")  # For demonstration, print a message to the console
    # open an executable file from the ./structure/ directory
    
    
    structureDirectory = os.path.join(os.path.dirname(__file__), 'structure')
    executableStructure = os.path.join(structureDirectory, 'structure.exe')  # Adjust the name as needed
    mainparamsFile = os.path.join(structureDirectory, 'mainparams')
    extraparamsFile = os.path.join(structureDirectory, 'extraparams')
    if os.path.exists(executableStructure):
        try:
            inputFile = os.path.join(structureDirectory, "input")
            with open(inputFile, 'w', encoding='utf-8') as input_file:
                input_file.write(structureFile)  # Write the opened structureFile content to input.txt
            
            # Change the working directory to ./structure/ and run the executable
            result = subprocess.run([executableStructure], cwd=structureDirectory, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

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


def editMainParameters():
    structureDirectory = os.path.join(os.path.dirname(__file__), 'structure')
    mainparamsFile = os.path.join(structureDirectory, 'mainparams')

    # Parameter labels and default values
    params = {
        "MAXPOPS (K; number of populations assumed)": "2",
        "BURNIN (length of burnin period)": "1000",
        "NUMREPS (MCMC reps after burnin)": "10000",
        "REPS (number of replicates)": "1",
        "PLOIDY": "2",
        "LABEL (1=yes, 0=no)": "1",
        "POPDATA (1=yes, 0=no)": "1",
        "POPFLAG (1=yes, 0=no)": "0",
        "LOCDATA (1=yes, 0=no)": "0",
        "PHENOTYPE (1=yes, 0=no)": "0",
        "EXTRACOLS (extra columns before genotype)": "0",
        "MARKERNAMES (1=yes, 0=no)": "1",
        "RECESSIVEALLELES (1=yes, 0=no)": "0",
        "MAPDISTANCES (1=yes, 0=no)": "0"
    }

    def submit():
        # Collect values
        user_values = {key: entry.get() for key, entry in entries.items()}
        try:
            with open(mainparamsFile, 'r', encoding='utf-8') as file:
                lines = file.readlines()

            new_lines = []
            for line in lines:
                updated = False
                for key in user_values:
                    define_key = key.split()[0]
                    if line.strip().startswith(f"#define {define_key}"):
                        new_line = f"#define {define_key} {user_values[key]}   // {key.split('(', 1)[-1]}\n"
                        new_lines.append(new_line)
                        updated = True
                        break
                if not updated:
                    new_lines.append(line)

            with open(mainparamsFile, 'w', encoding='utf-8') as file:
                file.writelines(new_lines)

            messagebox.showinfo("Success", "Parameters updated successfully.")
            window.destroy()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to update file: {e}")

    # Create window
    window = tk.Toplevel()
    window.title("Edit STRUCTURE Parameters")
    window.configure(bg='lightblue')

    tk.Label(window, text="Edit STRUCTURE Parameters", bg='lightblue',
             font=("Helvetica", 14, "bold")).pack(pady=10)

    form_frame = tk.Frame(window, bg='lightblue')
    form_frame.pack(padx=20, pady=10)

    entries = {}
    for i, (label, default) in enumerate(params.items()):
        tk.Label(form_frame, text=label, bg='lightblue', anchor='w').grid(row=i, column=0, sticky='w', pady=2)
        entry = tk.Entry(form_frame, width=25)
        entry.insert(0, default)
        entry.grid(row=i, column=1, pady=2)
        entries[label] = entry

    tk.Button(window, text="Submit", command=submit, bg='white').pack(pady=10)

def editExtraParameters():
    structureDirectory = os.path.join(os.path.dirname(__file__), 'structure')
    extraparamsFile = os.path.join(structureDirectory, 'extraparams')

    extra_params = {
    "NOADMIX (0=admixture, 1=no-admixture)": "0",
    "LINKAGE (use linkage model)": "0",
    "USEPOPINFO (use prior population info)": "0",
    "LOCPRIOR (use location info to improve weak data)": "0",
    "FREQSCORR (allele frequencies are correlated among pops)": "1",
    "ONEFST (same value of Fst for all)": "0",
    "INFERALPHA (infer admixture parameter)": "1",
    "POPALPHAS (individual alpha per population)": "0",
    "ALPHA (initial Dirichlet parameter)": "1.0",
    "INFERLAMBDA (infer lambda)": "0",
    "POPSPECIFICLAMBDA (separate lambda for each pop)": "0",
    "LAMBDA (allele frequency Dirichlet parameter)": "1.0"
    }

    def save_parameters(entries):
        try:
            with open(extraparamsFile, 'w', encoding='utf-8') as file:
                for label, entry in entries.items():
                    param_name = label.split()[0]
                    value = entry.get()
                    file.write(f"#define {param_name} {value} // {label}\n")
            result_label.config(text="Extra parameters saved successfully!", fg="green")
        except Exception as e:
            result_label.config(text=f"Error: {e}", fg="red")

    # GUI setup
    window = tk.Tk()
    window.title("Edit Extra Parameters - Program Options Only")
    window.configure(bg='light blue')

    frame = tk.Frame(window, bg='light blue')
    frame.pack(padx=20, pady=20)

    entries = {}
    for idx, (label, default) in enumerate(extra_params.items()):
        tk.Label(frame, text=label, bg='light blue', anchor='w', width=55).grid(row=idx, column=0, sticky='w')
        entry = tk.Entry(frame, width=10)
        entry.insert(0, default)
        entry.grid(row=idx, column=1)
        entries[label] = entry

    tk.Button(window, text="Save Parameters", command=lambda: save_parameters(entries)).pack(pady=10)

    result_label = tk.Label(window, text="", bg='light blue')
    result_label.pack()

    window.mainloop()   
    
