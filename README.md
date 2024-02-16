# Mortha Laboratory 
Codes used for image analysis 


#Code 1 (Figure number):
This code analyzes cellular interactions by calculating the distances between ILC2 or CD11b cells and B or T cells. It initially filters single-cell data for relevant phenotypes in both control and diseased conditions. Next, it measures the distance to the nearest cell from each phenotype in both conditions and combines this information with the original data. Subsequently, it generates a summary of the average distances to specific phenotypes for both conditions. Finally, it creates a density plot illustrating the distribution of distances

#Code 2 (Figure number):
This code conducts a comprehensive analysis of cellular interactions within healthy and infected conditions. Initially, it calculates the distances between key cell types, including T & B, B & ILC2, and ILC2 & T cells, generating tables summarizing these distances. Subsequently, the code identifies triplet interactions and filters duplicate occurrences. Through an iterative process, it evaluates the number of triplets where all three cell types are within a specified radius, providing insights into spatial relationships. Finally, the code visualizes these findings by plotting threshold values against the distance of T cells. Overall, the code serves to elucidate intricate cellular dynamics and spatial patterns in both control and diseased contexts.

#Code 3 (Figure number):
This code conducts an extensive analysis of cellular interactions under both healthy and infected conditions. Initially, it computes the distances between crucial cell types—T & B, B & ILC2, and ILC2 & T cells—and compiles tables summarizing these measurements. Subsequently, it identifies and filters duplicate triplet interactions. It then calculates the average triplet distance for each condition separately. 
