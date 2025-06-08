import tkinter as tk
import numpy as np
import pandas as pd
import structure2Genalex as s2g
import structureFunctions as sf
from collections import defaultdict, Counter 

################
# Target for genalex Functions
    # N – Number of individuals
    # Na – Number of alleles
    # Ne – Effective number of alleles
    # I – Shannon’s Information Index
    # Ho – Observed heterozygosity
    # He – Expected heterozygosity
    # UHe – Unbiased expected heterozygosity
    # FIS – Inbreeding coefficient
    # %P – Percentage of polymorphic loci

    # FST – Genetic differentiation among populations
    # GST – Nei’s genetic differentiation
    # Dest – Jost's D (allelic differentiation)
    # AMOVA – Analysis of Molecular Variance (with ΦST)
    # Pairwise FST / ΦST / Dest – Between each pair of populations

def populationStatistics():
    global namesOfLoci
    global numberOfLoci
    global numberOfIndividuals
    global numberOfPopulations
    global individualsPerPopn
    global structureDF
    global genalexDF


    print("\nComputing population statistics...\n")
    
##### PREPARING DATAFRAMES AND RESULTS TABLE
    # Create a DataFrame from the structure dataset
    computeDF = s2g.structureDF.copy()  # Create a copy of the structure DataFrame for GenAlEx format to not meddle the original structure DataFrame
    #print("\nStructure (computeDF) DataFrame: \n", computeDF)  # Print the computeDF for debugging
    # Create copy of genalex DataFrame for computations
    genalexCopy = s2g.genalexDF.copy()  # Create a copy of the GenAlEx DataFrame for computations
    #print("\nGenAlEx (genalexCopy) DataFrame: \n", genalexCopy)  # Print the copy for debugging
    # Reformat tab delimiters in genalexCopy to "/" 
    genalexSlash = genalexCopy.replace('\t', '/', regex=True)  # Replace tab delimiters with "/"
    #print("\nGenAlEx (slash-delimited) DataFrame: \n", genalexSlash)  # Print the copy after replacing tabs for debugging
    # Initialize a DataFrame to store results
    resultTable = pd.DataFrame(columns=["Population", "N", "Ne", "Na", "Ea" , "I", "Pa" , "Ho", "He", "Ht", "Fst", "Fis", "Fit"])
    populations = genalexCopy['Population'].unique()
    #print("\nResult Table initialized: \n", resultTable)  # Print the result table for debugging

##### RENAMING ALLELES FOR CONVENIENCE
    for locus in sf.namesOfLoci:
        if locus not in computeDF.columns:
            raise ValueError(f"Locus {locus} not found in the structure DataFrame. Please check the structure file.")
        unique_alleles = computeDF[locus].dropna().unique()
        # Exclude '-1' and '0' values from unique alleles
        unique_alleles = [allele for allele in unique_alleles if allele not in ['-1', '0']]  # Exclude '-1' and '0' values
        if len(unique_alleles) == 0:
            raise ValueError(f"No alleles found for locus {locus}. Please check the structure file.")
        # Sort the unique alleles by value, convert to integer, sort, and convert back to string
        unique_alleles = [int(allele) for allele in unique_alleles if allele.isdigit()]  # Convert to int and filter out non-digit alleles
        unique_alleles.sort()
        #print(f"\nUnique alleles for {locus} (Before renaming): {unique_alleles}")  # Print unique alleles for debugging
        # Rename the unique alleles to 1, 2, 3, etc.
        allele_mapping = {str(allele): str(i + 1) for i, allele in enumerate(unique_alleles)}
        #print(f"\nAllele mapping for {locus}: {allele_mapping}")  # Print the allele mapping for debugging
        # Apply the mapping to the structure DataFrame
        computeDF[locus] = computeDF[locus].map(allele_mapping)

##### COUNTING INDIVIDUALS PER POPULATION USING GENALLEX COPY
    print("\nCounting individuals per population ...")  # Debug message
    for population in populations:
        pop_data = genalexCopy[genalexCopy['Population'] == population]
        n = len(pop_data)  # Number of individuals
        if n == 0 :
            print(f"\nSkipping population {population} due to zero individuals or alleles.")
            continue  # Skip if no individuals or alleles
        resultTable.at[population, 'Population'] = population  # Add population name to result table
        resultTable.at[population, 'N'] = n  # Add number of individuals to result table
    print("\nResult Table after adding population and N: \n", resultTable)  # Print the result table for debugging
 
##### COMPUTE MEAN NUMBER OF ALLELES PER LOCUS PER POPULATION USING COMPUTEDF
    for population in computeDF['Population'].unique():
        pop_data = computeDF[computeDF['Population'] == population]
        if pop_data.empty:
            print(f"\nWarning: No data found for population {population}. Skipping this population.")
            continue    # Skip if no data for the population
        # Count unique alleles per locus
        na = pop_data.iloc[:, 2:].nunique()  # Count unique alleles per locus
        # Calculate mean number of alleles per locus
        mean_na = na.mean()
        print(f"\nMean number of alleles per locus for population {population}: {mean_na}")  # Print mean number of alleles for debugging
        # Add to result table
        resultTable.at[population, 'Na'] = mean_na  # Add mean number of alleles to result table
    print("\nResult Table after adding mean Na: \n", resultTable)  # Print the result table for debugging
        
    # Print the structure DataFrame after renaming alleles for debugging
    #print(f"\nStructure DataFrame after renaming alleles for {locus}: \n", computeDF)  # Print the structure DataFrame after renaming alleles for debugging

##### IDENTIFYING PRIVATE ALLELES USING GENALEX COPY (SLASH-DELIMITED genalexSlash)
    locusColumns = genalexSlash.columns[2:]
    # Count allele occurrences per population per locus
    allelePresence = defaultdict(lambda: defaultdict(set))

    for _, row in genalexSlash.iterrows():
        pop = str(row["Population"])  # Ensure consistent type
        for locus in locusColumns:
            genotype = row[locus]
            if isinstance(genotype, str) and "/" in genotype:
                allele1, allele2 = genotype.split("/")
                if allele1 != "0":
                    allelePresence[locus][allele1].add(pop)
                if allele2 != "0":
                    allelePresence[locus][allele2].add(pop)
    # Identify private alleles
    privateAlleles = defaultdict(lambda: defaultdict(list))  # private_alleles[pop][locus] = [alleles]

    for locus in locusColumns:
        for allele, pops in allelePresence[locus].items():
            allele = str(allele)  # Ensure allele is a string
            if len(pops) == 1 and allele != "0" and allele != "-1":
                pop = next(iter(pops))
                privateAlleles[pop][locus].append(allele)


    # Display results
    #from pprint import pprint
    #print("\nPrivate alleles per population and locus:")
    #pprint(privateAlleles)  # Print private alleles for debugging

    # Get count of private alleles per population
    print("Number of Loci SF: " + str(sf.numberOfLoci))
    for pop,loci in privateAlleles.items():
        # Get length of private alleles list for the population        
        for locus, alleles in loci.items():
            numberPrivateAlleles = len(alleles)  # Length of private alleles list for the population
        #   print(f"\nPopulation {pop} has {numberPrivateAlleles} private alleles for locus {locus}.")  # Print private allele count for debugging
        # Compute mean number of private alleles per population with total number of loci as denominator
        
        meanPrivateAlleles = (np.sum([len(alleles) for alleles in loci.values()]))/sf.numberOfLoci  # Mean number of private alleles per population
        #print(f"\nPopulation {pop} has {meanPrivateAlleles} mean private alleles across all loci.")  # Print mean private allele count for debugging
        # Add private allele counts to result table
        resultTable.at[pop, 'Pa'] = meanPrivateAlleles
    print("\nResult Table after adding mean Pa (private alleles): \n", resultTable)  # Print the result table for debugging  

##### Make empty dictionary to store Shannon's Information Index (I) values per population
    shannonIndex = {}  # Dictionary to store Shannon's Information Index (I) values per population
    # print("\nInitialized empty dictionary for Shannon's Information Index (I): ", shannonIndex)  # Print the initialized dictionary for debugging

##### COMPUTE ALLELE FREQUENCIES, EXP HETEROZYGOSITY, and EFFECTIVE NUMBER OF ALLELES AND WRITING TO FILE
    # Create a dataframe for locus within populations
    expHeterozygosityDF = pd.DataFrame(index=populations, columns=computeDF.columns[2:])  # Initialize empty DataFrame for expected heterozygosity
    #print("\nInitialized empty DataFrame for expected heterozygosity: \n", expHeterozygosityDF)  # Print the initialized DataFrame for debugging

    # Initialize a dictionary to store dataframes for allele frequencies
    alleleFreqDF = {}  # Dictionary to store allele frequencies DataFrames for each locus

    with open("allele.freq", "w") as f:
        for locus in computeDF.columns[2:]:
            print(f"\nProcessing locus: {locus}")  # Optional debug

            # Count allele occurrences per population
            alleleCounts = computeDF.groupby('Population')[locus].value_counts().unstack(fill_value=0)
            # print(f"\nAllele counts for {locus}:\n", alleleCounts)  # Print allele counts for debugging
            if alleleCounts.empty:
                print(f"\nWarning: No allele counts found for locus {locus}. Skipping this locus.")
                continue

            # Sort allele columns numerically
            sorted_columns = sorted(alleleCounts.columns, key=lambda x: int(x))
            alleleCounts = alleleCounts[sorted_columns]

            # print(f"\nSorted allele counts for {locus}:\n", alleleCounts)  # Print sorted allele counts for debugging
            # Compute frequencies and round
            alleleFrequencies = alleleCounts.div(alleleCounts.sum(axis=1), axis=0).round(4)
            #print(f"\nAllele frequencies for {locus}:\n", alleleFrequencies)  # Print allele frequencies for debugging
            squaredAlleleFrequencies = alleleFrequencies ** 2  # Square the allele frequencies
            # Compute 1 - sum of squared frequencies per population
            #print(f"\nSquared allele frequencies for {locus}:\n", squaredAlleleFrequencies)  # Print squared allele frequencies for debugging
            expHeterozygosity = 1 - squaredAlleleFrequencies.sum(axis=1)  # Expected heterozygosity per population

            for population in expHeterozygosity.index:
                expHeterozygosityDF.at[population, locus] = expHeterozygosity.loc[population]
            #print(f"\nExpected heterozygosity for {locus}:\n", expHeterozygosity)

            # Write locus label
            locus_number = list(computeDF.columns[2:]).index(locus) + 1  # nth number of the locus (1-based)
            f.write(f"@Locus {locus_number}: {locus}\n")
            #print(f"\nAllele frequencies for {locus}:\n", alleleFrequencies)  # Print allele frequencies for debugging

##### Compute Shannon's Information Index (I)
            for population in alleleFrequencies.index:
                p = alleleFrequencies.loc[population]
                i = -np.nansum(p * np.log(p + 1e-9))
                #print(f"\nShannon's Information Index (I) for {population} at {locus}: {i}")
                # store in dictionary
                shannonIndex[(population, locus)] = i       

            # Write frequencies with tab separation, no extra blank lines
            alleleFrequencies.to_csv(f, sep='\t', lineterminator='\n')
            f.write("\n")
            # Store alleleFrequencies in allelefreqDF dictionary
            alleleFreqDF[locus] = alleleFrequencies
            
    print("\nAllele frequencies written to allele.freq file.")  # Confirmation message
    # print("\nShannon's Information Index (I) values per population: ", shannonIndex)  # Print the Shannon's Information Index for debugging

    print("\nAllele Frequencies DataFrame: \n", alleleFreqDF)  # Print the allele frequencies DataFrame for debugging


    # Get mean Shannon's Information Index (I) for each population
    for population in populations:
        # Filter the dictionary for the current population
        pop_indices = [key for key in shannonIndex.keys() if key[0] == population]
        if not pop_indices:
            print(f"\nWarning: No Shannon's Information Index (I) found for population {population}.")
            continue  # Skip if no indices found for the population
        # Calculate mean I for the population
        mean_i = np.mean([shannonIndex[key] for key in pop_indices])
        resultTable.at[population, 'I'] = mean_i  # Add mean I to result table

    print("\nResult Table after adding Shannon's Information Index (I): \n", resultTable)  # Print the result table for debugging
    #print("\nExpected heterozygosity DataFrame: \n", expHeterozygosityDF)  # Print the expected heterozygosity DataFrame for debugging

##### Compute mean expected heterozygosity (He) for each population in expHeterozygosityDF
    for population in populations:
        if population not in expHeterozygosityDF.index:
            print(f"\nWarning: Population {population} not found in expHeterozygosityDF. Skipping this population.")
            continue
        meanHe = expHeterozygosityDF.loc[population].mean()
        #print(f"\nMean expected heterozygosity (He) for population {population}: {meanHe}")
        resultTable.at[population, 'He'] = meanHe  # Add mean He to result table
    print("\nResult Table after adding expected heterozygosity (He): \n", resultTable)  # Print the result table for debugging

##### Calculate effective number of alleles (Ne) for each population using effNeDF and expHeterozygosityDF
    effNeDF = pd.DataFrame(index=populations, columns=computeDF.columns[2:])  # Initialize empty DataFrame for effective number of alleles
    #print("\nInitialized empty DataFrame for effective number of alleles: \n", effNeDF)  # Print the initialized DataFrame for debugging
    



##### COMPUTING HETEROZYGOSITY USING GENALEX COPY (SLASH-DELIMITED genalexSlash)
    print("\nComputing heterozygosity and homozygosity ...\n")  # Debug message
    #print(genalexSlash)  # Print columns of genalexSlash for debugging
    # Count genotypes 
    heterozygous = defaultdict(lambda: defaultdict(int))  # heterozygosity[population][locus] = count
    homozygous = defaultdict(lambda: defaultdict(int))  # homozygosity[population][locus] = count

    for _, row in genalexSlash.iterrows():
        population = row['Population']
        for locus in locusColumns:
            genotype = row[locus]
            if isinstance(genotype, str) and "/" in genotype:
                alleles = genotype.split("/")
                if len(alleles) == 2 and all(a not in ("0", "-1") for a in alleles):  # Ensure valid genotype
                    # Count heterozygous and homozygous genotypes
                    if alleles[0] != alleles[1]:
                        heterozygous[population][locus] += 1
                    else:
                        homozygous[population][locus] += 1

##### Create a dataframe to store heterozygous for each population by locus
    heterozygousDF = pd.DataFrame(index=populations, columns=locusColumns)
    # Store heterozygous[population][locus] counts in the heterozygousDF DataFrame
    for population in populations:
        for locus in locusColumns:
            heterozygousDF.at[population, locus] = heterozygous[population][locus]
    #print("\nHeterozygous DataFrame: \n", heterozygousDF)  # Print the heterozygous DataFrame for debugging

    # Compute observed heterozygosity (Ho) using heterozygousDF per population per locus
    heterozygosityDF = pd.DataFrame(index=populations, columns=locusColumns)
    #print("\nInitialized empty heterozygosity DataFrame: \n", heterozygosityDF)  # Print the initialized DataFrame for debugging

    # Refer to resultTable for individuals per population
    for population in populations:
        if population not in resultTable.index:
            print(f"\nWarning: Population {population} not found in resultTable. Skipping this population.")
            continue
        for locus in locusColumns:
            if locus not in heterozygousDF.columns:
                print(f"\nWarning: Locus {locus} not found in heterozygousDF. Skipping this locus.")
                continue
            # Calculate observed heterozygosity
            observedHeterozygosity = heterozygousDF.at[population, locus] / resultTable.at[population, 'N']  # Divide heterozygous count by total count for the population
            #print(f"\nObserved heterozygosity for population {population} at locus {locus}: {observedHeterozygosity}")  # Print observed heterozygosity for debugging
            if pd.isna(observedHeterozygosity):
                print(f"\nWarning: No observed heterozygosity found for population {population} at locus {locus}. Skipping this locus.")
                continue
            heterozygosityDF.at[population, locus] = observedHeterozygosity  # Store observed heterozygosity in the DataFrame
    #print("\nHeterozygosity DataFrame: \n", heterozygosityDF)  # Print the heterozygosity DataFrame for debugging

    # Calculate mean observed heterozygosity (Ho) for each population
    for population in populations:
        if population not in heterozygosityDF.index:
            print(f"\nWarning: Population {population} not found in heterozygosityDF. Skipping this population.")
            continue
        meanHo = heterozygosityDF.loc[population].mean()  # Calculate mean observed heterozygosity for the population
        #print(f"\nMean observed heterozygosity (Ho) for population {population}: {meanHo}")  # Print mean Ho for debugging
        resultTable.at[population, 'Ho'] = meanHo  # Add mean Ho to result table
    print("\nResult Table after adding observed heterozygosity (Ho): \n", resultTable)  # Print the result table for debugging


    
##### COMPUTING F-STATISTICS
    # Create the allele frequency DataFrame

    # # Compute Fst
    # for population in populations:
    #     if population not in heterozygosityDF.index:
    #         print(f"\nWarning: Population {population} not found in heterozygosityDF. Skipping this population.")
    #         continue
        
    #     # Subtract He in population from Ht to get Fst
    #     print(resultTable.at[population, "He"])  # Print He for debugging

    #     Fst = (Ht- resultTable.at[population, "He"]) / Ht  # Fst calculation 

    #     print(f"\nFst for population {population}: {Fst}")
    #     # add Fst to result table
    #     resultTable.at[population, 'Fst'] = Fst
    #     Fis = (resultTable.at[population, 'He'] - resultTable.at[population, 'Ho']) / resultTable.at[population, 'He']
    #     resultTable.at[population, 'Fis'] = Fis
    #     Fit = (resultTable.at[population, 'Ht'] - resultTable.at[population, 'Ho']) / resultTable.at[population, 'Ht']
    #     resultTable.at[population, 'Fit'] = Fit

    # print("\nResult Table after computing Fst, Fis, & Fit: \n", resultTable)  # Print the result table for debugging


    #  genalexSlash

    #     # Calculate effective number of alleles (Ne)
    #     allele_counts = pop_data.iloc[:, 2:].apply(lambda x: x.value_counts(), axis=0).fillna(0)
    #     ne = (1 / (allele_counts.apply(lambda x: (x / x.sum()) ** 2).sum(axis=0)).sum() if na > 0 else 0)

    #     # Calculate Shannon's Information Index (I)
    #     p = allele_counts.div(allele_counts.sum(axis=1), axis=0)
    #     i = -(p * np.log(p + 1e-9)).sum(axis=1).mean() if not p.empty else 0

    #     # Calculate observed heterozygosity (Ho)
    #     ho = pop_data.apply(lambda row: sum(row[2:].str.split('\t').apply(lambda x: len(set(x)) > 1)), axis=1).mean()

    #     # Calculate expected heterozygosity (He)
    #     he = 2 * (allele_counts / allele_counts.sum(axis=1, keepdims=True)).apply(lambda x: x ** 2).sum(axis=1).mean()

    #     results = results.append({
    #         "Population": population,
    #         "N": n,
    #         "Na": na,
    #         "Ne": ne,
    #         "I": i,
    #         "Ho": ho,
    #         "He": he
    #     }, ignore_index=True)

        




        




    