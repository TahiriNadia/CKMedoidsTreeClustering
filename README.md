# Clustering_Phylogenetic_trees
A new fast method for building multiple consensus trees using k-medoids

# About
	=> =====================================================================================
	=> Program : K-Medoids super-trees - 2015
	=> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)
	=> This program computes a clustering for phylogenetic trees based on the k-Medoids.
	=> =====================================================================================

# Installation
	$ make 
	$ make install
	$ make clean
	$ make help

# Examples
	Please execute next command line:
	=> For simulation: ./cSuperTree -simulation number_species percent_nose number_of_simulation int_parameter
	=> For trees: ./cSuperTree -tree nameFile_trees int_parameter
	=> For matrice: ./cSuperTree -matrice nameFile_matrice int_parameter

	List of criterion with K-medoid algorithm
	=> option 1 - Calinski-Harabasz with RF
	=> option 2 - Calinski-Harabasz with RF squared
	=> option 3 - Silhouette with RF
	=> option 4 - Silhouette with RF squared

	where RF is Robinson and Foulds distance
	
	./cSuperTree -tree data_47.txt 1
	
# Input
	Phylogenetic trees file in Netwick format (see example data/data_47.txt)
	
# Output
	stat.csv : statitic values (score of CH or SH), partitionning data and K found
	output.txt : each phylogenetic tree was attribued of their cluster 
