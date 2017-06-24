//================================================================================================
//=  HGT-DETECTION v3.3.1
//=  Modified by : Nadia Tahiri
//=  Date of modification : March 2014
//=  Authors : Alix Boc and Vladimir Makarenkov
//=  Date : November 2009
//=  
//=  Description : This program detect horizontal gene transfer (HGT). As input it takes 2 
//=  trees: a species tree and a gene tree. the goal is to transform the species tree
//=  into the gene tree following a transfer scenario. There are 3 criteria : the robinson and
//=  Foulds distance, the least-square criterion and the bipartition distance. We also use the
//=  subtree constraint. With this version we can now perform simulation.
//=
//=	 input   : file with species tree and gene tree in the newick format.
//=			   In case of simulation, the species tree and all the gene trees in the same file in
//=            the phylip format or newick string
//=  output  : a list of HGT and the criteria values for each one.
//=	 options :
//=
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#include <iostream>
using namespace std;
#pragma warning(disable:4996)

#include "structures.h"
#include "utils_tree.cpp"
#include "fonctions.cpp"

#define binaireSpecies 0 
#define binaireGene    1
int compteur=0;


void traiterSignal(int sig){
	printf("\nMESSAGE : SEGMENTATION FAULT #%d DETECTED",sig);
	printf("\nUse valgrind or gdb to fix the problem");
	printf("\n");
	exit(-1);
}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
void main_hgt(string tree1, string tree2, double *distances ){

	struct InputTree SpeciesTree;				    //== initial species tree
	/* struct InputTree SpeciesTreeCurrent;		//== initial species tree */
/* 	struct InputTree FirstTree; */
	/* struct InputTree FirstTree2; */
	/* struct InputTree FirstTreeB; */
	struct InputTree GeneTree;					    //== initial gene tree
	
	struct InputTree SpeciesTreeReduce;				    //== initial species tree reduit
	struct InputTree GeneTreeReduce;					    //== initial gene tree reduit
	
	struct InputTree SpeciesTreeRed;			  //== reduced species tree
	struct InputTree GeneTreeRed;				    //== reduced gene tree
	/* struct ReduceTrace aMap;					      //== mapping structure between species tree and recuded species tree */
	/* struct InputTree geneTreeSave; */
	/* struct HGT * bestHGTRed = NULL;				  //== list of HGT for the reduced tree */
	/* struct HGT * bestHGT = NULL;				    //== list of HGT for the normale tree */
	/* struct HGT * outHGT = NULL; */
	/* struct HGT * bestHGTmulticheck = NULL; */
	/* int nbHGT_boot; */
	/* int first = 1,k,l; */
	int cpt_hgt,i,j,nb_same_espece,nbTree=0;
	int min_diff = 0; // difference minimum of species between T1 and nb_same_espece or between T2 and nb_same_espece
	int bootstrap = 0;
	int multigene = 0;
	int nbHgtFound = 0;
	/* struct CRITERIA * multicheckTab=NULL; */
	struct CRITERIA aCrit;						       //== struture of all the criteria
/* 	struct DescTree *DTSpecies,					       //== structure of submatrices for the species tree
					*DTGene;	 */				       //== structure of submatrices for the gene tree
	struct Parameters param;
	/* FILE *in,*out; */
	/* int max_hgt,nbHGT; */
	/* int ktSpecies; */
    	/* int trivial = 1; */
	
	/* int *speciesLeaves = NULL; */
	/* int RFref; */
	/* int imc;	 */
	/* char *mot = (char*)malloc(100); */
	/* int nomorehgt=0; */
	/* initInputTree(&geneTreeSave); */

	//== read parameters
	if(readParameters(&param)==-1){
		printf("\nhgt : no options specified, see the README file for more details\n");
		exit(-1);
	}

	rand_bootstrap = param.rand_bootstrap;
	
	if(strcmp(param.speciesroot,"file") == 0){
		if(!file_exists(param.speciesRootfileLeaves) && !file_exists(param.speciesRootfile)){
			printf("\nhgt : The file %s does not exist",param.speciesRootfileLeaves);
			exit(-1);
		}
	}
	if(strcmp(param.generoot,"file") == 0){
		if(!file_exists(param.geneRootfileLeaves) && !file_exists(param.geneRootfile)){
			printf("\nhgt : The file %s does not exist",param.geneRootfileLeaves);
			exit(-1);
		}
	}
		
/* 	initInputTree(&FirstTree);
	initInputTree(&SpeciesTreeCurrent); */
	
//==============================================================================
//============================= LECTURE DES ARBRES =============================
//==============================================================================
	initInputTree(&SpeciesTree);
	initInputTree(&GeneTree);
	nb_same_espece = readInputFile(tree1,tree2, param.input,&SpeciesTree,&GeneTree,param.errorFile);
	
	if(nb_same_espece<0) {
		nb_same_espece=0;
	}
/* 	if(nb_same_espece==-1) {
		printf("\nCannot read input data !!\n");
		exit(-1);
	} */
  
/* 	cpt_hgt = 0; */

/* 	initInputTree(&SpeciesTree);
	initInputTree(&GeneTree);
	initInputTree(&SpeciesTreeRed);
	initInputTree(&GeneTreeRed); */

	initInputTree(&SpeciesTreeReduce);
	initInputTree(&GeneTreeReduce);
	//== lecture des matrices ou chaines newick en entree
	if(readInput(SPECIE,param.input,&SpeciesTreeReduce) == -1){ printf("\nError in species tree\n"); exit(-1);}
	if(readInput(GENE,param.input,&GeneTreeReduce) == -1){ printf("\nError in gene tree\n"); getchar(); exit(-1);}

	TrierMatrices(GeneTreeReduce.Input,GeneTreeReduce.SpeciesName,SpeciesTreeReduce.SpeciesName,SpeciesTreeReduce.size);
	
	NJ(SpeciesTreeReduce.Input,SpeciesTreeReduce.ADD,SpeciesTreeReduce.size);
	NJ(GeneTreeReduce.Input,GeneTreeReduce.ADD,GeneTreeReduce.size);

	//== construction des differentes repr�sentation des arbres (adjacence,aretes,longueur,degre)
	CreateSubStructures(&SpeciesTreeReduce,1,binaireSpecies);
	CreateSubStructures(&GeneTreeReduce,1,binaireGene);

	InitCriteria(&aCrit,SpeciesTreeReduce.size);
	computeCriteria(SpeciesTreeReduce.ADD,GeneTreeReduce.ADD,SpeciesTreeReduce.size,&aCrit,SpeciesTreeReduce.LONGUEUR,SpeciesTreeReduce.ARETE,GeneTreeReduce.LONGUEUR,GeneTreeReduce.ARETE);
	distances[0]=aCrit.RF;
	if(SpeciesTree.size-nb_same_espece>GeneTree.size-nb_same_espece){
		min_diff = GeneTree.size-nb_same_espece;
	}else{
		min_diff = SpeciesTree.size-nb_same_espece;
	}
	min_diff = fabs(min_diff);
	min_diff = min_diff * min_diff;
	
/* 	cout<<"SpeciesTree "<<SpeciesTree.size<<endl;
	cout<<"GeneTree "<<GeneTree.size<<endl;
	cout<<"nb_same_espece "<<nb_same_espece<<endl;
	cout<<"MIN DIFF "<<min_diff<<endl;
	cout<<"distances[0] "<<distances[0]<<endl; */
	
	if(nb_same_espece==3){
		distances[0]=0;
	}else{
		/* distances[0]=(distances[0]/(2.0*nb_same_espece-6.0))*((SpeciesTree.size+GeneTree.size+min_diff)/(SpeciesTree.size+GeneTree.size)); //normalisation par le nombre d'espèces communes entre les deux arbres */
		distances[0]=distances[0]/(2.0*nb_same_espece-6.0); //normalisation par le nombre d'espèces communes entre les deux arbres
	}
	/* cout<<distances[0]<<endl; */
	distances[1]=aCrit.LS;
	distances[2]=aCrit.BD;
	distances[3]=nb_same_espece;	
	
/* 	for(int i=0;i<2*SpeciesTreeReduce.size;i++){
		free(SpeciesTreeReduce.ADD[i]);
		free(SpeciesTreeReduce.Adjacence[i]);
		free(SpeciesTreeReduce.Input[i]);
		free(SpeciesTreeReduce.W[i]);
	}
	
	for(int i=0;i<=SpeciesTreeReduce.size;i++){
		free(SpeciesTreeReduce.SpeciesName[i]);
	} */
	
/* 	free(SpeciesTreeReduce.degre);
	free(SpeciesTreeReduce.ADD);
	free(SpeciesTreeReduce.Adjacence);
	free(SpeciesTreeReduce.Input);
	free(SpeciesTreeReduce.W);
	free(SpeciesTreeReduce.ARETE);
	free(SpeciesTreeReduce.LONGUEUR);
	free(SpeciesTreeReduce.SpeciesName); */

	
	FreeCriteria(&aCrit,SpeciesTreeReduce.size);
	freeInputTree(&SpeciesTree,SpeciesTree.size);
	freeInputTree(&GeneTree,GeneTree.size);
	freeReducedTree(&SpeciesTreeReduce,SpeciesTreeReduce.size);
	freeReducedTree(&GeneTreeReduce,GeneTreeReduce.size);
	
/* 	freeReducedTree(&SpeciesTree,SpeciesTree.size);
	freeReducedTree(&GeneTree,GeneTree.size); */
	
/* 	cout<<"taille finale "<<nb_same_espece<<endl;
 	cout<<"\nCOMPTEUR == "<<compteur++<<endl;
	cout<<"distance RF "<<distances[0]<<endl;
	cout<<"distance LS "<<distances[1]<<endl;
	cout<<"distance BD "<<distances[2]<<endl; */
	/* return distances; */
}
