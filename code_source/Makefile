all: KMedoidsTreeClustering.cpp
	g++ KMedoidsTreeClustering.cpp -openmp -g -o KMTC
	chmod -R 777 * .
	@echo "Please execute the following command line:"
	@echo '===> For trees: ./KMTC -tree Input_file int_parameter'
	@echo ' '
	@echo '===> ================================================================================================================================================='
	@echo '===> Program   : KMedoidsTreeClustering - 2017'
	@echo '===> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)'
	@echo '===> This program computes a clustering for phylogenetic trees based on the K-Medoid partitioning.'
	@echo '===> This program computes a clustering of phylogenetic trees based on the K-Medoids partitioning algorithm.'
	@echo '===> The set of trees in the Newick format should be provided as input.'
	@echo '===> The optimal partitioning in K classes is returned as output. The number of classes can be determined by the '
	@echo '===> Silhouette and Calinski–Harabasz cluster validity indices adapted for tree clustering. The non-squared and '
	@echo '===> squared Robinson and Foulds topological distance can be used. The recommended option: Silhouette + non-squared Robinson and Foulds distance.'
	@echo '===> ================================================================================================================================================='
	@echo ' '
	
install: KMedoidsTreeClustering.cpp
	g++ KMedoidsTreeClustering.cpp -openmp -g -o KMTC
	chmod -R 777 * .
	@echo "Please execute the following command line:"
	@echo '===> For trees: ./KMTC -tree Input_file int_parameter'
	@echo ' '
	@echo '===> ================================================================================================================================================='
	@echo '===> Program   : KMedoidsTreeClustering - 2017'
	@echo '===> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)'
	@echo '===> This program computes a clustering for phylogenetic trees based on the K-Medoid partitioning.'
	@echo '===> This program computes a clustering of phylogenetic trees based on the K-Medoids partitioning algorithm.'
	@echo '===> The set of trees in the Newick format should be provided as input.'
	@echo '===> The optimal partitioning in K classes is returned as output. The number of classes can be determined by the '
	@echo '===> Silhouette and Calinski–Harabasz cluster validity indices adapted for tree clustering. The non-squared and '
	@echo '===> squared Robinson and Foulds topological distance can be used. The recommended option: Silhouette + non-squared Robinson and Foulds distance.'
	@echo '===> ================================================================================================================================================='
	@echo ' '
	
clean: 
	$(RM) KMTC
	
help: 
	@echo 'Please execute the following command lines (on two steps):'
	@echo '===> 1) make install'
	@echo '===> 2) For trees: ./KMTC -tree Input_file int_parameter'
	@echo ' '
	@echo '===> ================================================================================================================================================='
	@echo '===> Program   : KMedoidsTreeClustering - 2017'
	@echo '===> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)'
	@echo '===> This program computes a clustering for phylogenetic trees based on the K-Medoid partitioning.'
	@echo '===> This program computes a clustering of phylogenetic trees based on the K-Medoids partitioning algorithm.'
	@echo '===> The set of trees in the Newick format should be provided as input.'
	@echo '===> The optimal partitioning in K classes is returned as output. The number of classes can be determined by the '
	@echo '===> Silhouette and Calinski–Harabasz cluster validity indices adapted for tree clustering. The non-squared and '
	@echo '===> squared Robinson and Foulds topological distance can be used. The recommended option: Silhouette + non-squared Robinson and Foulds distance.'
	@echo '===> ================================================================================================================================================='
	@echo ' '
	@echo 'List of criteria with K-medoid algorithm'
	@echo 'option 1 - Calinski-Harabasz with RF'
	@echo 'option 2 - Calinski-Harabasz with RF squared'
	@echo 'option 3 - Silhouette with RF'
	@echo 'option 4 - Silhouette with RF squared'
