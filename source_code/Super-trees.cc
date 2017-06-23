#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#include "K-medoids.cpp"

using namespace std;

int main_consense(char **argv,vector<int> tabIndices, vector <string> mesTrees,int intParam){
	 // Variables
    time_t tbegin,tend;
    double texec=0.;

    // Start timer
    tbegin=time(NULL);				// get the current calendar time
	
/* 	printf ("Consensus'tree K-means partitioning\n");
	printf("Nadia Tahiri and Vladimir Makarenkov - Departement d'informatique - Universite du Quebec a Montreal\n");
	printf ("Original code by :  Pierre Legendre - Departement de sciences biologiques - Universite de Montreal.\n");
	printf ("(c) Pierre Legendre, 1999\n");
 */
	//Varriables
	double **Matrice_RF;
	double **Ww;
	double **n_identique;
	
	double *distances = new double[4];
	string tree1;
	string tree2;
	
	for (int j=0; j<4; j++){
		distances[j]=0.0;
	}
				
	//Création de la matrice carrée et symétrique (mesTrees.size()*mesTrees.size()) : Matrice_RF
	Matrice_RF= new double*[mesTrees.size()];
	Ww= new double*[mesTrees.size()];
	n_identique= new double*[mesTrees.size()];
	
	
	for(int lineDist=0;lineDist<mesTrees.size();lineDist++){
		Matrice_RF[lineDist]= new double[mesTrees.size()];
		Ww[lineDist]= new double[mesTrees.size()];
		n_identique[lineDist]= new double[mesTrees.size()];
	}
	
	
	double MoinsUn = -1.0;
	double PlusUn = 1.0;
	double Tirage1 = 0.0;
	double Tirage2 = 0.0;
	
	// Remplissage de la matrice des distances RF en faisant appel à main_hgt 
	// qui calcule la distance RF entre chaque paire d'arbre
	for(int line=0;line<mesTrees.size();line++){
		//mettre des valeurs 0.0 pour la diagonale de la matrice RF
		Matrice_RF[line][line]=0.0;
		
		for(int column=line;column<mesTrees.size();column++){

			// Affectation de deux arbres : tree1 et tree2
			tree1=mesTrees[line];
			tree2=mesTrees[column];
			// Appel des algorithmes des calcules des distances : RF
			main_hgt(tree1,tree2,distances);
			Matrice_RF[line][column]=distances[0];
			
			if(Matrice_RF[line][column]<=0.0){
				Matrice_RF[line][column]=0.0;
			}
			
			if(isnan(Matrice_RF[line][column])){
				Matrice_RF[line][column]=0.0;
			}

			//Recuperer le nombre d'espèces communes
			n_identique[line][column]=distances[3];
			//Et la symétrie
			n_identique[column][line]=distances[3];	
			
			// pour remplir la symétrique de la matrice RF sans réaliser de calcul (car matrice carrée symétrique)
			Matrice_RF[column][line]=Matrice_RF[line][column];		

		}
	}
	
	//creation de la matrice de distances RF : mat
	 for (int i=0; i<mesTrees.size(); i++)
	 {
		 for (int j=0; j<mesTrees.size(); j++)
		 {
			Ww[i][j]=1.0;
		 }
	 }
	 
	//appel de l'algorithme de K-means:
		if(mesTrees.size()>7){	
			main_kmedoids(argv,mesTrees,Matrice_RF,n_identique,Ww,tabIndices,intParam);
		}	

	
	//Liberation of memory
	for (int i=0;i<mesTrees.size();i++){
		delete [] Matrice_RF[i];
		delete [] Ww[i];
		delete [] n_identique[i];
	}
	delete [] Matrice_RF;
	delete [] Ww;
	delete [] n_identique;
	delete [] distances;
	// End timer
    tend=time(NULL);				// get the current calendar time

	// Compute execution time
    texec=difftime(tend,tbegin);	// tend-tbegin (result in second)
	return 0;
}
