#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hgt_int.cpp"
#include "Super-trees.cc"


using namespace std;

struct TNode{	
    int NoNoeud;
	char *seq;
	TNode * gauche;
	TNode * droit;
	int nbFils;
};


//==Prototypes des fonctions
void SAVEASNewick(double *LONGUEUR, long int *ARETE,int nn,const char* t); 
void createTree(int nbSpecies,double l,const char* t);
void GenererMatrice(TNode **TabNoeud, double **DISS);
void Tree_edges(double **,long int *,double *);
int floor1(double x);
void odp1ct(double**,int*,int*,int*,int);
void tree_generation(double **DA, double **DI, int n, double Sigma);
string longueur (string tree2);
int ExtraireDonneesLC(const char * chaine, char *champs, char * contenu);
void validation(int &intParam);


 #define INFINI 999999.99
 #define MaxRF 0
 double MaxLong=0;
 double MinLong = INFINI;
 double seuil;
 double epsilona = 0.00005;
 int n;
 
 void presenterProgramme(){
	printf ("KMTC - Clustering phylogenetic trees using k-medoids\n");
	printf("Nadia Tahiri and Vladimir Makarenkov - Departement d'informatique - Universite du Quebec a Montreal\n");
	printf("Choose the criterion to perform the k-medoids clustering of trees:\n");
	printf("option 1 - Calinski-Harabasz with RF\n");
	printf("option 2 - Calinski-Harabasz with RF squared\n");
	printf("option 3 - Silhouette with RF\n");
	printf("option 4 - Silhouette with RF squared\n");
 }


//=============================================================
//=========================== MAIN ============================
//=============================================================
int main(int nargs,char ** argv){

	presenterProgramme();
	printf("Traitement in progress... \n");
	
	char champs[100];
	char contenu[100];
		
	if(nargs < 2){
		printf("\nbad input..\nusage:%s {-tree}\n",argv[0]);
		exit(1);
	}		
	
	if(ExtraireDonneesLC(argv[1],champs,contenu)==1){
		if(strcmp("tree",champs) == 0){
			if(nargs != 4){
				printf("\nbad input..\nusage:%s {-tree} nameFile Parametre\n",argv[0]);
				exit(1);
			}	
			fstream fichier(argv[2]);
			int intParam = atoi(argv[3]);
			validation(intParam);
			vector <string> mesTrees;
			int ligne = 1;
			char ** cl2 = new char*[4];
			for (int i=0;i<4;i++){
				cl2[i] = new char[10];
			}
						
			strcpy(cl2[0], "*");
			strcpy(cl2[1], "?");
			strcpy(cl2[2], "?");
			strcpy(cl2[3], "?");
			vector <int> tabIndices;
			
			//read input file (tree file)
			if( !fichier ){
				cout << "File "<<argv[2]<<" no exist."<<endl;
			}else{
				while( !fichier.eof()){
					mesTrees.push_back("");//creation d'une ligne vide
					getline(fichier, mesTrees.back());//lecture d'une ligne du fichier
					ligne = mesTrees.size() - 1;//je recupere la taille du tableau (-1 pour la ligne 0)
					
					if(mesTrees[ligne].empty())//si la ligne est vide
						mesTrees.pop_back();//on la retire du tableau
				} 
				tabIndices.push_back(mesTrees.size());
				
				//call the method which compute Robinson and Foulds distance between each pair of trees
				main_consense(cl2,tabIndices,mesTrees,intParam);
				
				//vider les vecteurs
				mesTrees.clear();
				tabIndices.clear();
			}
		}		
	}
	printf("end of computation\n");
	return 0;
}

string longueur (string tree2){
	//en cours........		
	string lg ="";
	int cp = 0;
	int pours = 0;
	
	while(cp<tree2.length() && pours==0){
		if(tree2.at(cp)==':'){
			cp++;
			while((tree2.at(cp)>='0' && tree2.at(cp)<='9') || tree2.at(cp)>='.'){
				lg += tree2.at(cp);
				cp++;
			}
			pours = 1;
		}
		cp++;
	}	
	
	return lg;
}

//=============================================
//
//=============================================
void createTree (int nbSpecies,double l, const char * outputfilename){	 
    int i,j;
    double *LONGUEUR,Sigma;        //= Longueurs des aretes
    long int *ARETE;               //= Liste des aretes de l'arbre
    double **DA,**DI;		       //= Distance d'arbre
    
    n = nbSpecies;
    
	// Allocation dynamique pour les matrices
	DI = (double **) malloc((n+1)*sizeof(double*));
	DA=(double **) malloc((2*n-1)*sizeof(double*));

	for (i=0;i<=n;i++){
		DI[i]=(double*)malloc((n+1)*sizeof(double));
		 
		if (DI[i]==NULL)
		{
			printf(" Data matrix is too large(2)\n"); 
			exit(1);
		}
	}

    for(i=0;i<=2*n-2;i++){
         DA[i]=(double*)malloc((2*n-1)*sizeof(double));
    }	 

    ARETE=(long int *) malloc((4*n)*sizeof(long int));
    LONGUEUR=(double *) malloc((2*n)*sizeof(double));

	//===== CA COMMENCE ICI =====
	
	
	//==================================
	//= CREATION D'UN ARBRE ALEATOIRE
	//==================================
	tree_generation(DA, DI, n, Sigma = 0.0);
	Tree_edges(DA, ARETE, LONGUEUR);
	
	MaxLong = 0;
	MinLong = INFINI;
    
	for(int j=1;j<=2*n-3;j++){
		MaxLong += LONGUEUR[j-1];
	} 			
    //printf("\nlongueur moyenne des branches = %lf",MaxLong/(2*n-3));
    
    double facteur = 1.0/(MaxLong/(2*n-3));
    MaxLong = 0;
    
    for(int j=1;j<=2*n-3;j++){
        LONGUEUR[j-1] = (LONGUEUR[j-1]*facteur)*l;
		MaxLong += LONGUEUR[j-1];
	} 			
	//printf("\nl= %lf,facteur = %lf,longueur moyenne des branches = %lf\n",l,facteur,MaxLong/(2*n-3));
    //=============================================================================
	//= SAUVEGARDE DE L'ARBRE
	//=============================================================================
    SAVEASNewick(LONGUEUR,ARETE,n,outputfilename); 

	
	for (i=0;i<=n;i++){
		free(DI[i]);
	}
	free(DI);
    for(i=0;i<=2*n-2;i++){
         free(DA[i]);
    }	
	free(DA);
	free(ARETE);
	free(LONGUEUR);

	
}

//=================================
//
//=================================
double P(double l){
	return 0.25*(1.0+3.0*exp(-4.0*1*l)); 
}

// This C++ function is meant to generate a random tree distance matrix DA
// of size (nxn) and a distance (i.e. dissimilarity) matrix DI obtained
// from DA by adding a random normally distributed noise with mean 0 
// and standard deviation Sigma. For more detail, see Makarenkov and
// Legendre (2004), Journal of Computational Biology.  
void tree_generation(double **DA, double **DI, int n, double Sigma)
{
   struct TABLEAU { int V; } **NUM, **A;
   int i,j,k,p,a,a1,a2,*L,*L1,n1;
   double *LON,X0,X,U;
   //time_t t;
   
   n1=n*(n-1)/2;
    
   L=(int *)malloc((2*n-2)*sizeof(int));
   L1=(int *)malloc((2*n-2)*sizeof(int));
   LON=(double *)malloc((2*n-2)*sizeof(double));
    
   NUM=(TABLEAU **)malloc((2*n-2)*sizeof(TABLEAU*));
   A=(TABLEAU **)malloc((n1+1)*sizeof(TABLEAU*));
    
   for (i=0;i<=n1;i++)
   {
      A[i]=(TABLEAU*)malloc((2*n-2)*sizeof(TABLEAU));
      if (i<=2*n-3) NUM[i]=(TABLEAU*)malloc((n+1)*sizeof(TABLEAU));
      
      if ((A[i]==NULL)||((i<=2*n-3)&&(NUM[i]==NULL))) 
      {
       printf("\nData matrix is too large\n");
       exit(1);
      }
    }  

/* Generation of a random additive tree topology T*/

   for (j=1;j<=2*n-3;j++)
   {
    for (i=1;i<=n;i++)
     {
       A[i][j].V=0;
       NUM[j][i].V=0;
     }
    for (i=n+1;i<=n1;i++)
      A[i][j].V=0;
   }
   A[1][1].V=1; L[1]=1; L1[1]=2; NUM[1][1].V=1; NUM[1][2].V=0;
     
//   srand((unsigned) time(&t));

   for (k=2;k<=n-1;k++)
   {
    p=(rand() % (2*k-3))+1;
    for (i=1;i<=(n*(k-2)-(k-1)*(k-2)/2+1);i++)
     A[i][2*k-2].V=A[i][p].V;
    for (i=1;i<=k;i++)
    {
     a=n*(i-1)-i*(i-1)/2+k+1-i;
     if (NUM[p][i].V==0)
      A[a][2*k-2].V=1;
     else
      A[a][p].V=1;
    }
    for (i=1;i<=k;i++)
    {
      a=n*(i-1)-i*(i-1)/2+k+1-i;
      A[a][2*k-1].V=1;
    }
    for (j=1;j<=k;j++)
    {
     if (j==L[p])
     {
       for (i=1;i<=2*k-3;i++)
       {
        if (i!=p)
        {
         if (L1[p]>L[p])
         a=floor1((n-0.5*L[p])*(L[p]-1)+L1[p]-L[p]);
         else
         a=floor1((n-0.5*L1[p])*(L1[p]-1)+L[p]-L1[p]);
         if (A[a][i].V==1)
         {
          if (NUM[i][L[p]].V==0)
            a=floor1((n-0.5*L[p])*(L[p]-1)+k+1-L[p]);
          else
            a=floor1((n-0.5*L1[p])*(L1[p]-1)+k+1-L1[p]);
          A[a][i].V=1;
         }
        }
       }
     }
     else if (j!=L1[p])
     {
      a=floor1((n-0.5*j)*(j-1)+k+1-j);
      if (j<L[p])
      a1=floor1((n-0.5*j)*(j-1)+L[p]-j);
      else
      a1=floor1((n-0.5*L[p])*(L[p]-1)+j-L[p]);

      if (j<L1[p])
      a2=floor1((n-0.5*j)*(j-1)+L1[p]-j);
      else
      a2=floor1((n-0.5*L1[p])*(L1[p]-1)+j-L1[p]);
      for (i=1;i<=2*k-3;i++)
       {
        if ((i!=p)&&((A[a1][i].V+A[a2][i].V==2)||((NUM[i][j].V+NUM[i][L[p]].V==0)&&(A[a2][i].V==1))||((NUM[i][j].V+NUM[i][L1[p]].V==0)&&(A[a1][i].V==1))))
          A[a][i].V=1;
       }
     }
    }
    for (i=1;i<=k;i++)
     NUM[2*k-2][i].V=NUM[p][i].V;
    NUM[2*k-2][k+1].V=1;

    for (i=1;i<=k;i++)
     NUM[2*k-1][i].V=1;

    for (i=1;i<=2*k-3;i++)
    {
     if (((NUM[i][L[p]].V+NUM[i][L1[p]].V)!=0)&&(i!=p))
      NUM[i][k+1].V=1;
    }

    L[2*k-2]=k+1; L1[2*k-2]=L1[p];
    L[2*k-1]=L1[p]; L1[2*k-1]=k+1;
    L1[p]=k+1;
   }

   //srand((unsigned) time(&t));
   
 // Exponential distribution generator with expectation 1/(2n-3)  
 // 1. Generate U~U(0,1)
 // 2. Compute X = -ln(U)/lambda
 // 3. To obtain an exponential distribution with theoretical mean (= "expected
 // value") of 1/(2n-3), use lambda = 1 and multiply X by 1/(2n-3):     
 // Thus, X = -(1/(2n-3))*ln(U)
 // This formula is given in Numerical Recipes, p.200.*/ 
 // Every branch length of the tree T is then multiplied by 1+aX,
 // where X followed the standard exponential distribution (P(X>n)=exp(-n)) and 
 //  "a" is a tuning factor to adjust the deviation from the molecular clock; "a" is set at 0.8. 
 // 4. LON[i] = LON[i]*(1+0.8*-ln(U))
 
   U = 0.1;
   while (U<0.1) U = 1.0*rand()/RAND_MAX;
   
   for (i=1;i<=2*n-3;i++)
		LON[i] = -1.0/(2*n-3)*log(U);    
   
   //for (i=1;i<=2*n-3;i++)
   i=1;
   while (i<=2*n-3)
   {  U = 1.0*rand()/RAND_MAX;  LON[i] = 1.0*LON[i]*(1.0+0.8*(-log(U)));  LON[i] = LON[i]*0.8; ///5.0 ;//0.8
           if (LON[i]>2*epsilona) i++;/*printf("U=%f LON[%d]=%f \n",U,i,LON[i]);*/}
      
 // Computation of a tree distance matrix (tree metric matrix)
   for (i=1;i<=n;i++)
   {
    DA[i][i]=0;
    for (j=i+1;j<=n;j++)
    {
     DA[i][j]=0;
     a=floor1((n-0.5*i)*(i-1)+j-i);
     for (k=1;k<=2*n-3;k++)
      if (A[a][k].V==1) DA[i][j]=DA[i][j]+LON[k];
     DA[j][i]=DA[i][j];
    }
   }

// Normalisation of the tree distance matrix

/* Addition of noise to the tree metric */
  //srand((unsigned) time(&t));
  
  for (i=1;i<=n;i++)
  {
   DI[i][i]=0.0;
   for (j=i+1;j<=n;j++)
   {
    X=0.0;   
    for (k=1;k<=5;k++)
    {
	  X0 = 1.0*rand()/RAND_MAX;
      X=X+0.0001*X0;
    }
    X=2*sqrt(0.6)*(X-2.5);
    U=X-0.01*(3*X-X*X*X);
    DI[i][j]=DA[i][j]+Sigma*U;
    if (DI[i][j]<0) DI[i][j]=0.01;
    DI[j][i]=DI[i][j];
   }
  }
   
  free (L);
  free(L1);
  free(LON);    
  for (i=0;i<=n1;i++)
  {
    free(A[i]);
    if (i<=2*n-3) free(NUM[i]);
  }
  free(NUM);
  free(A);
}


int floor1(double x)
{  
  int i;
  
  if (ceil(x)-floor(x)==2) i=(int)x; 
  else if (fabs(x-floor(x)) > fabs(x-ceil(x))) i=(int)ceil(x);
  else i=(int)floor(x);
  return i;
} 

/****************************************************
* application de la m�thode NJ pour construire
* une matrice additive
****************************************************/
void NJ(double **D1,double **DA)
{
double **D,*T1,*S,*LP,Som,Smin,Sij,L,Lii,Ljj,l1,l2,l3,epsilona = 0.00005;
int *T,i,j,ii,jj,n1;

D=(double **) malloc((n+1)*sizeof(double*));
T1=(double *) malloc((n+1)*sizeof(double));
S=(double *) malloc((n+1)*sizeof(double));
LP=(double *) malloc((n+1)*sizeof(double));
T=(int *) malloc((n+1)*sizeof(int));

for (i=0;i<=n;i++)
{
  D[i]=(double*)malloc((n+1)*sizeof(double));

  if (D[i]==NULL)
  {
	  { printf("Data matrix is too large"); return;}
  }
}


L=0;
Som=0;
for (i=1;i<=n;i++)
{
  S[i]=0; LP[i]=0;
  for (j=1;j<=n;j++)
    {
      D[i][j]=D1[i][j];
      S[i]=S[i]+D[i][j];
    }
  Som=Som+S[i]/2;
  T[i]=i;
  T1[i]=0;
}

/* Procedure principale */
for (n1=n;n1>3;n1--)
{

  /* Recherche des plus proches voisins */
  Smin=2*Som;
  for (i=1;i<=n1-1;i++)
  {
   for (j=i+1;j<=n1;j++)
   {
    Sij=2*Som-S[i]-S[j]+D[i][j]*(n1-2);
    if (Sij<Smin)
    {
     Smin=Sij;
     ii=i;
     jj=j;
    }
   }
  }
/* Nouveau groupement */

  Lii=(D[ii][jj]+(S[ii]-S[jj])/(n1-2))/2-LP[ii];
  Ljj=(D[ii][jj]+(S[jj]-S[ii])/(n1-2))/2-LP[jj];

/* Mise a jour de D */

  if (Lii<2*epsilona) Lii=2*epsilona;
  if (Ljj<2*epsilona) Ljj=2*epsilona;
  L=L+Lii+Ljj;
  LP[ii]=0.5*D[ii][jj];

  Som=Som-(S[ii]+S[jj])/2;
  for (i=1;i<=n1;i++)
  {
    if ((i!=ii)&&(i!=jj))
    {
     S[i]=S[i]-0.5*(D[i][ii]+D[i][jj]);
     D[i][ii]=(D[i][ii]+D[i][jj])/2;
     D[ii][i]=D[i][ii];
    }
  }
  D[ii][ii]=0;
  S[ii]=0.5*(S[ii]+S[jj])-D[ii][jj];

  if (jj!=n1)
  {
    for (i=1;i<=n1-1;i++)
    {
     D[i][jj]=D[i][n1];
     D[jj][i]=D[n1][i];
    }
    D[jj][jj]=0;
    S[jj]=S[n1];
    LP[jj]=LP[n1];
  }
/* Mise a jour de DA */
  for (i=1;i<=n;i++)
  {
    if (T[i]==ii) T1[i]=T1[i]+Lii;
    if (T[i]==jj) T1[i]=T1[i]+Ljj;
  }


  for (j=1;j<=n;j++)
  {
   if (T[j]==jj)
   {
     for (i=1;i<=n;i++)
     {
       if (T[i]==ii)
       {
         DA[i][j]=T1[i]+T1[j];
         DA[j][i]=DA[i][j];
       }
     }
   }
  }

  for (j=1;j<=n;j++)
  if (T[j]==jj)  T[j]=ii;

  if (jj!=n1)
  {
   for (j=1;j<=n;j++)
   if (T[j]==n1) T[j]=jj;
  }
}

/*Il reste 3 sommets */

l1=(D[1][2]+D[1][3]-D[2][3])/2-LP[1];
l2=(D[1][2]+D[2][3]-D[1][3])/2-LP[2];
l3=(D[1][3]+D[2][3]-D[1][2])/2-LP[3];
if (l1<2*epsilona) l1=2*epsilona;
if (l2<2*epsilona) l2=2*epsilona;
if (l3<2*epsilona) l3=2*epsilona;
L=L+l1+l2+l3;

for (j=1;j<=n;j++)
{
   for (i=1;i<=n;i++)
   {
    if ((T[j]==1)&&(T[i]==2))
    {
     DA[i][j]=T1[i]+T1[j]+l1+l2;
     DA[j][i]=DA[i][j];
    }
    if ((T[j]==1)&&(T[i]==3))
    {
     DA[i][j]=T1[i]+T1[j]+l1+l3;
     DA[j][i]=DA[i][j];
    }
    if ((T[j]==2)&&(T[i]==3))
    {
     DA[i][j]=T1[i]+T1[j]+l2+l3;
     DA[j][i]=DA[i][j];
    }
   }
   DA[j][j]=0;
}

  free(T);
  free(T1);
  free(S);
  free(LP);
  for (i=0;i<=n;i++)
  {
  	free(D[i]);
  }
  free(D);

}


/* La recherche d'un ordre diagonal plan */
void odp(double **D,int *X,int *i1,int *j1)
{
 double S1,S;
 int i,j,k,a,*Y1;
 
 Y1=(int *) malloc((n+1)*sizeof(int));
 
 for(i=1;i<=n;i++)
  Y1[i]=1;
  
 X[1]=*i1;
 X[n]=*j1;
 if (n==2) return;
 Y1[*i1]=0;
 Y1[*j1]=0;
 for(i=0;i<=n-3;i++)
 { a=2;
   S=0;
   for(j=1;j<=n;j++)
   { if (Y1[j]>0)
    {
      S1= D[X[n-i]][X[1]]-D[j][X[1]]+D[X[n-i]][j];
      if ((a==2)||(S1<=S))
      {
        S=S1;
        a=1;
        X[n-i-1]=j;
        k=j;
      }
    }
   }
   Y1[k]=0;
 }
 
   free(Y1);
}

/* Calcule les tableaux ARETE et LONGUEUR contenant les aretes et leur longueurs
a partir d'une matrice de distance d'arbre DI de taille n par n; 
la fonction auxiliaire odp est appell�e une fois */  

void Tree_edges (double **DI, long int *ARETE, double *LONGUEUR)
{ 
 struct EDGE { unsigned int U; unsigned int V; double LN;};
 struct EDGE *Path,*Tree;
 int i,j,k,p,P,*X;
 double S,DIS,DIS1,*L,**D;
 
  X=(int *)malloc((n+1)*sizeof(int));  
  L=(double *)malloc((n+1)*sizeof(double));
  Tree=(EDGE *)malloc((2*n-2)*sizeof(EDGE));
  Path=(EDGE *)malloc((n+2)*sizeof(EDGE));
 
 
 D=(double **) malloc((n+1)*sizeof(double*));
 
 for (i=0;i<=n;i++)
 {
  D[i]=(double*)malloc((n+1)*sizeof(double)); 
  
  if (D[i]==NULL)
  {
    printf("Data matrix is too large"); exit(1); 
  }
 }
 
 i=1; j=n;
 odp1ct(DI,X,&i,&j,n);
 
 for (i=1;i<=n;i++)
 { 
  for (j=1;j<=n;j++) 
   D[i][j]=DI[i][j];
 }  
   
 /* Verification de l'equivalence des topologies */
 L[1]=D[X[1]][X[2]];
 Path[1].U=X[1];
 Path[1].V=X[2];

 p=0;
 P=1;
 
 for(k=2;k<=n-1;k++)
 {

  DIS=(D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
  DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
  
  S=0.0;
  i=0;

  if (DIS>2*epsilona)
  {
    while (S<DIS-epsilona)
    {
     i=i+1;
     S=S+L[i];
    }
  }
  else { DIS=0; i=1; }
  
  Tree[p+1].U=n+k-1;
  Tree[p+1].V=Path[i].V;
  Tree[p+1].LN=S-DIS;
  if (Tree[p+1].LN<2*epsilona) Tree[p+1].LN=epsilona; 
  
  for (j=i+1;j<=P;j++)
  {
   Tree[p+j-i+1].U=Path[j].U;
   Tree[p+j-i+1].V=Path[j].V;
   Tree[p+j-i+1].LN=L[j];
   if (L[j]<2*epsilona) L[j]=2*epsilona; 
  }
  p=p+P-i+1;
  
  Path[i].V=n+k-1;
  Path[i+1].U=n+k-1;
  Path[i+1].V=X[k+1];
  L[i]=L[i]+DIS-S;
  L[i+1]=DIS1;
  P=i+1;
 }
 
  for (i=1;i<=P;i++) 
  {
   Tree[p+i].U=Path[i].U;
   Tree[p+i].V=Path[i].V;
   Tree[p+i].LN=L[i];
  }
 
 for (i=1;i<=2*n-3;i++)
 {
  if (fabs(Tree[i].LN-epsilona)<=2*epsilona)
   Tree[i].LN=0.0;
  ARETE[2*i-2]=Tree[i].U;
  ARETE[2*i-1]=Tree[i].V;
  LONGUEUR[i-1]=Tree[i].LN;   
  if (LONGUEUR[i-1]<2*epsilona) LONGUEUR[i-1] = 2*epsilona;
 } 
 
 free(X);
 free(Tree);
 free(L);
 free(Path);

 
 for (i=0;i<=n;i++)    
    free(D[i]);
    
 free(D);
}

/* Fonction auxiliaire odp. La recherche d'un ordre diagonal plan */

void odp1ct(double **D, int *X, int *i1, int *j1, int n)
{
 double S1,S;
 int i,j,k,a,*Y1;
 
 Y1=(int *) malloc((n+1)*sizeof(int));
 
 for(i=1;i<=n;i++)
  Y1[i]=1;
  
 X[1]=*i1;
 X[n]=*j1;
 if (n==2) return;
 Y1[*i1]=0;
 Y1[*j1]=0;
 for(i=0;i<=n-3;i++)
 { a=2;
   S=0;
   for(j=1;j<=n;j++)
   { if (Y1[j]>0)
    {
      S1= D[X[n-i]][X[1]]-D[j][X[1]]+D[X[n-i]][j];
      if ((a==2)||(S1<=S))
      {
        S=S1;
        a=1;
        X[n-i-1]=j;
        k=j;
      }
    }
   }
   Y1[k]=0;
 }     
   free(Y1);
}


void SAVEASNewick(double *LONGUEUR, long int *ARETE,int nn,const char* t) 
{
	int n,root,a;
	int Ns;
	int i, j, sd, sf, *Suc, *Fre, *Tree, *degre, *Mark;
	double *Long;
	int *boot; 
	char *string = (char*)malloc(100000);	
	n = nn;
	Ns=2*n-3;
	
	double * bootStrap= NULL;
	
	Suc =(int*) malloc((2*n) * sizeof(int));
	Fre =(int*) malloc((2*n) * sizeof(int));
	degre =(int*) malloc((2*n) * sizeof(int));
	Long = (double*) malloc((2*n) * sizeof(double));	
	boot = (int*) malloc((2*n) * sizeof(int));
	Tree = (int*) malloc((2*n) * sizeof(int));
	Mark =(int*) malloc((2*n) * sizeof(int));
	
	if ((degre==NULL)||(Mark==NULL)||(string==NULL)||(Suc==NULL)||(Fre==NULL)||(Long==NULL)||(Tree==NULL)||(ARETE==NULL)||(LONGUEUR==NULL))	
		{ printf("Tree is too large to be saved"); return;} 
	
	for (i=1;i<=2*n-3;i++)
	{ 
		
		if (i<=n) degre[i]=1;
		else degre[i]=3;
	} degre[2*n-2]=3;
	
	root=Ns+1;
	
	int cpt=0;
	
	for (;;)
	{
		cpt++;
		if(cpt > 1000000) exit(1);
		a=0; a++;
		for (j=1;j<=2*n-2;j++)
			Mark[j]=0;
		
		for (i=1;i<=2*n-3;i++)
		{ 	  									
			if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-2]] = (int) bootStrap[i-1];
				
			}
			else if ((degre[ARETE[2*i-1]]==1)&&(degre[ARETE[2*i-2]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-1]]=ARETE[2*i-2]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-1]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-1]] = (int) bootStrap[i-1];
				
			}
			else if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]==1)&&(Mark[ARETE[2*i-2]]==0)&&(Mark[ARETE[2*i-1]]==0))
			{
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; root=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; a=-1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-2]] = (int) bootStrap[i-1];
			}
			if (a==-1) break;
		}
		if (a==-1) break;
	}
	
	
	/*  On decale et on complete la structure d'arbre avec Successeurs et Freres  */
	for (i=Ns+1;i>0;i--)
	{ 	Fre[i]=0; Suc[i]=0;
	}	Tree[root]=0;/*Tree[Ns+1]=0;*/
	
	for (i=1;i<=Ns+1/*Ns*/;i++)
	{	
		if (i!=root) 
		{
			sd=i; sf=Tree[i];
			if (Suc[sf]==0) Suc[sf]=sd;
			else {	
				sf=Suc[sf];
				while (Fre[sf]>0) sf=Fre[sf];
				Fre[sf]=sd;
			}		 
		}
	}
	
	
	/* On compose la chaine parenthesee */
	string[0]=0; i=root;/*i=Ns+1;*/
	cpt=0;
	for (;;)
	{	
		if(cpt > 1000000) exit(1);
		
		if (Suc[i]>0)
		{	sprintf(string,"%s(",string);
		Suc[i]=-Suc[i]; i=-Suc[i]; }
		else if (Fre[i]!=0)
		{	if (Suc[i]==0) sprintf(string,"%s%d:%.4f,",string,i,Long[i]);
			else {
				if(bootStrap != NULL)
					sprintf(string,"%s%d:%.4f,",string,boot[i],Long[i]);
				else
					sprintf(string,"%s:%.4f,",string,Long[i]);
			}
		i=Fre[i]; }
		else if (Tree[i]!=0)
		{	if (Suc[i]==0) sprintf(string,"%s%d:%.4f)",string,i,Long[i]);
		else {
			if(bootStrap != NULL)
				sprintf(string,"%s%d:%.4f)",string,boot[i],Long[i]);
			else
				sprintf(string,"%s:%.4f)",string,Long[i]);
			}
		i=Tree[i]; }
		else break;
	}	
	strcat(string,";");
	
	FILE *pt_t = fopen(t,"a+");
	fprintf(pt_t,"%s\n",string);
	fclose(pt_t);
	
	free(Suc); free(Fre); free(Tree); free(Long); free(degre); free(Mark);	free(string);
}

//=====================================================================
//
//=====================================================================
int ExtraireDonneesLC(const char * chaine, char *champs, char * contenu){

	int cpt=0;
	int tailleChaine;

	if(chaine[0] != '-'){
		return 0;
	}else{
		tailleChaine = (int)strlen(chaine);
		for(int i=0;i<tailleChaine;i++){
			champs[i] = chaine[i+1];
		}
	}

	return 1;
}


void validation(int &intParam){

	while(intParam<0 || intParam>19){
		cout<<endl;
		cout<<"Invalid Parameter. Please enter again. "<<endl;
		cout<<"1  : K-medoid (with criterion CH) "<<endl;
		cout<<"2  : K-medoid (with criterion CH and RF squared) "<<endl;
		cout<<"3  : K-medoid (with criterion SH) "<<endl;
		cout<<"4  : K-medoid (with criterion SH and RF squared) "<<endl;
		cout<<"0  : EXIT"<<endl;
		
		cin>>intParam;
	}
	
	if(intParam==0){
		cout<<"====END OF PROGRAM!===="<<endl;
		exit(0);
	}
}
