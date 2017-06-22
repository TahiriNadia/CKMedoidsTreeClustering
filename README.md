# Clustering_Phylogenetic_trees
A new fast method for building multiple consensus trees using k-medoids

# About
	=> =============================================================================================================
	=> Program   : KMedoidsTreeClustering - 2017
	=> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)
	=> This program computes a clustering of phylogenetic trees based on the K-Medoids partitioning algorithm.
	=> The set of trees in the Newick format should be provided as input.
	=> The optimal partitioning in K classes is returned as output. The number of classes can be determined by the 
	=> Silhouette and CalinskiâHarabasz cluster validity indices adapted for tree clustering. The non-squared and 
	=> squared Robinson and Foulds topological distance can be used. 
	=> The recommended option: Silhouette + non-squared Robinson and Foulds distance.
	=> =============================================================================================================

# Installation
	$ make 
	$ make install
	$ make clean
	$ make help

# Examples
	Please execute the following command line:
	=> For trees: ./KMTC -tree Input_file int_parameter

	List of criteria with K-medoid algorithm
	=> option 1 - Calinski-Harabasz with RF
	=> option 2 - Calinski-Harabasz with RF squared
	=> option 3 - Silhouette with RF
	=> option 4 - Silhouette with RF squared

	where RF is Robinson and Foulds distance
	
	./KMTC -tree ../data/trees/input.tre 3
	
# Input
	=> See the folder data
	Phylogenetic trees file in Netwick format (see example data/input.tre)
	
# Output
	=> See the folder Output
	1) stat.csv : statitic values (score of CH or SH), partitionning data and K found
	2) output.txt : each phylogenetic tree was attribued of their cluster 
