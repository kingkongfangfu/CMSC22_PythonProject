import structureFunctions as sf
import pandas as pd
import numpy as np


def struc2Genalex():
    """
    Convert a structure file to GenAlEx format.
    
    Parameters:
    structureFile (str): Path to the input structure file.
    outputFile (str): Path to the output GenAlEx file.
    
    Returns:
    None
    """
    global structureFile
    global structureDataset
    global namesOfLoci
    global numberOfLoci
    global numberOfIndividuals
    global numberOfPopulations
    global individualsPerPopn
    global genalexDF
    global structureDF

    # Check if the structure file is provided
    if not sf.structureFile:
        raise ValueError("No structure file provided. Please select a file.")

    # Read the structure file
    #print(sf.structureFile)

    individualsPerLocus = "\t".join(str(len(value)) for value in sf.individualsPerPopn.values())
    print("\nIndividuals per locus: ", individualsPerLocus)  # Print the individuals per locus for debugging

    genalexHeader = str(sf.numberOfLoci) + "\t" + str(sf.numberOfIndividuals) + "\t" + str(sf.numberOfPopulations) + "\t" + individualsPerLocus
    genalexSubheader = "\nPopulation\tIndividual" + sf.structureFile[0] + "\t".join(sf.namesOfLoci)

    #print("\nHeader: ", genalexHeader)  # Print the header for debugging
    #print("\nSubheader: ", genalexSubheader)  # Print the subheader for debugging

    # Extract the relevant data from the structure file
    data = []
    listOfIndividuals = []

    for line in sf.structureDataset:
        parts = line.strip().split()
        
        # Check line format
        if len(parts) != sf.numberOfLoci + 2:
            print(f"\nWarning: Line does not match expected format: {line.strip()}")
            continue

        individual = parts[0]
        population = parts[1]
        alleles = parts[2:]

        if len(alleles) != sf.numberOfLoci:
            print(f"\nWarning: Number of alleles does not match expected number for individual {individual}: {line.strip()}")
            continue

        # Construct data row: [population, individual, allele1, allele2, ...]
        entry = [population, individual] + alleles

        if individual in listOfIndividuals:
            # Get index of existing entry and merge allele data
            index = listOfIndividuals.index(individual)
            for i in range(sf.numberOfLoci):
                data[index][i + 2] += "\t" + alleles[i]
        else:
            # First time seeing this individual
            listOfIndividuals.append(individual)
            data.append(entry)

    # Print the data for debugging
    #print("\nData: ", data)  # Print the data for debugging
    #print("\nList of individuals: ", listOfIndividuals)  # Print the list of individuals for debugging

    # Define the output GenAlEx file path
    genalexFile = "output_genalex_file.txt"

    # Write the data to a GenAlEx file
    with open(genalexFile, 'w') as f:
        f.write(genalexHeader)
        f.write("\n")  # Add a newline after the header
        f.write(genalexSubheader+ "\n")
        print("\nNumber of rows: ", len(data))  # Print the number of rows for debugging
        # Write the data rows
        for row in data:
            f.write("\t".join(row) + "\n")

        print("\nData written to GenAlEx file: ", genalexFile)  # Print the output file name for debugging

    # Create a DataFrame from the data
    print("\nCreating Genalex DataFrame from data...")  # Print a message for debugging
    genalexDF = pd.DataFrame(data, columns=["Population", "Individual"] + sf.namesOfLoci)
    print("\nGenAlEx DataFrame: \n", genalexDF)  # Print the DataFrame for debugging

    # Load or assume your DataFrame as `df`, with first two columns being 'Population' and 'Individual'
    locus_columns = genalexDF.columns[2:]  # All locus names
    hetero_matrix = pd.DataFrame(index=locus_columns, columns=sorted(genalexDF['Population'].unique()))

    for locus in locus_columns:
        for population in genalexDF['Population'].unique():
            subset = genalexDF[genalexDF['Population'] == population][locus].dropna()
            heterozygous = 0
            total = 0

            for genotype in subset:
                alleles = genotype.split('\t')
                if len(alleles) == 2:
                    total += 1
                    if alleles[0] != alleles[1]:
                        heterozygous += 1

            # Avoid division by zero
            freq = heterozygous / total if total > 0 else 0
            hetero_matrix.loc[locus, population] = round(freq, 4)

    # Optional: save to file
    hetero_matrix.to_csv("heterozygosity_per_locus.txt", sep="\t")


    # Making a dataframe from the structure dataset
    structureData = []
    for entry in sf.structureDataset:
        parts = entry.split()
        individual = parts[0]
        population = parts[1]
        alleles = parts[2:]
        structureData.append([population, individual] + alleles)

    # Create a DataFrame from the structure data
    structureDF = pd.DataFrame(structureData, columns=["Population", "Individual"] + sf.namesOfLoci)

    #print("\nStructure DataFrame: \n", structureDF)  # Print the DataFrame for debugging




