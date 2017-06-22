all: createTree_ConsensusTrees_simulations.cpp
	g++ createTree_ConsensusTrees_simulations.cpp -openmp -g -o cSuperTree
	chmod -R 777 * .
	@echo "Please execute next command line:"
	@echo '===> For trees: ./cSuperTree -tree nameFile_trees int_parameter'
	@echo ' '
	@echo '===> ========================================================================================================'
	@echo '===> Program : K-Medoid super-trees - 2017'
	@echo '===> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)'
	@echo '===> This program computes a clustering for phylogenetic trees based on the K-Medoid partitioning.'
	@echo '===> ========================================================================================================'
	@echo ' '
	
install: createTree_ConsensusTrees_simulations.cpp
	g++ createTree_ConsensusTrees_simulations.cpp -openmp -g -o cSuperTree
	chmod -R 777 * .
	@echo "Please execute next command line:"
	@echo '===> For trees: ./cSuperTree -tree nameFile_trees int_parameter'
	@echo ' '
	@echo '===> ========================================================================================================'
	@echo '===> Program : K-Medoid super-trees - 2017'
	@echo '===> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)'
	@echo '===> This program computes a clustering for phylogenetic trees based on the K-Medoid partitioning.'
	@echo '===> ========================================================================================================'
	@echo ' '
	
clean: 
	$(RM) cSuperTree
	
help: 
	@echo ' '
	@echo '===> ========================================================================================================'
	@echo '===> Program : K-Medoid super-trees - 2017'
	@echo '===> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)'
	@echo '===> This program computes a clustering for phylogenetic trees based on the K-Medoid partitioning.'
	@echo '===> ========================================================================================================'
	@echo ' '
	@echo "Please execute next command line:"
	@echo '===> For trees: ./cSuperTree -tree nameFile_trees int_parameter'
	@echo ' '
	@echo 'List of criteria with K-medoid algorithm'
	@echo 'option 1 - Calinski-Harabasz with RF'
	@echo 'option 2 - Calinski-Harabasz with RF squared'
	@echo 'option 3 - Silhouette with RF'
	@echo 'option 4 - Silhouette with RF squared'