//===================================================================================
//============================ DEFINITION DES FONCTIONS =============================
//===================================================================================

void deleteFiles(){
	
}

void trier(struct HGT *bestHGTRed, int nbHgtFound){
	struct HGT tmp;
	int inversion;
	int i;
	do{
		inversion=0;

		for(i=0;i<nbHgtFound-1;i++){
			if (bestHGTRed[i].crit.RF>bestHGTRed[i+1].crit.RF){
				tmp = bestHGTRed[i];
				bestHGTRed[i] = bestHGTRed[i+1];
				bestHGTRed[i+1] = tmp;
				inversion=1;
			}
		}
	}while(inversion);
}

//=================================================================
//== 
//=================================================================
void addLeafAndUpdate(struct InputTree *aTree, int choix){
	
	int n = aTree->size;
	int i,j;

	for(i=1;i<=2*n-3-aTree->kt;i++){
		if(aTree->ARETE[2*i-1] > n) aTree->ARETE[2*i-1] += 1 ;
		if(aTree->ARETE[2*i-2] > n) aTree->ARETE[2*i-2] += 1 ;
	}
	

	//= ajout des 2 nouvelles branches
	aTree->ARETE[2*(2*n-2-aTree->kt)-1] = aTree->ARETE[2*choix-1];
	aTree->ARETE[2*(2*n-2-aTree->kt)-2] = 2*n-aTree->kt;
	aTree->LONGUEUR[(2*n-2-aTree->kt)-1] = ((aTree->LONGUEUR[choix-1]/2.0)>5*epsilon)?(aTree->LONGUEUR[choix-1]/2.0):5*epsilon;

	aTree->ARETE[2*(2*n-1-aTree->kt)-1] = aTree->ARETE[2*choix-2];
	aTree->ARETE[2*(2*n-1-aTree->kt)-2] = 2*n-aTree->kt;
	aTree->LONGUEUR[(2*n-1-aTree->kt)-1] = ((aTree->LONGUEUR[choix-1]/2.0)>5*epsilon)?(aTree->LONGUEUR[choix-1]/2.0):5*epsilon;
	
	aTree->ARETE[2*choix-1] = n+1;
	aTree->ARETE[2*choix-2] = 2*n-aTree->kt;
	aTree->LONGUEUR[choix-1]= 5*epsilon;  

	aTree->size = aTree->size+1;

	loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,n+1,aTree->kt);

	Floyd(aTree->Adjacence,aTree->ADD,aTree->Input,n+1,aTree->kt);//5eme fois

	strcpy(aTree->SpeciesName[aTree->size],"Root");

}

void ListeSommets_taille_0(double ** matrix,int * tab_sommet,int size){
	
	int cpt=1;
	int cpt_init;
	int nouveau_cas=0;
	
	//== recherche du premier element admissible
	while(cpt <= size){
		if(tab_sommet[cpt] == 0)
			break;
		else
			cpt++;
	}
	
	if(cpt < size){
		//== recherche des distances nulles
		cpt_init = cpt++;
		while(cpt <= size){
			if(matrix[cpt_init][cpt] < epsilon ){
				tab_sommet[cpt] = 1;
				nouveau_cas = 1;
			}		
			cpt++;
		}
		if(nouveau_cas == 1){
			tab_sommet[cpt_init] = 1;
		}else
			tab_sommet[cpt_init] = 2;
	}else
		tab_sommet[0] = TRUE; //== aucun elements trouves
}

int Est_un_sous_ensemble_exact(int * sous_ensemble,int * ensemble){
	
	int temoin;
	
	if(sous_ensemble[0] > ensemble[0]){
		return 0;
	}
	else{
		temoin=0;
		for(int i=1;i<=sous_ensemble[0];i++){
			for(int j=1;j<=ensemble[0];j++){
				if(sous_ensemble[i] == ensemble[j]){
					temoin++;
				}
			}
		}
		if(temoin == sous_ensemble[0]) return 1;
	}
	return 0;
}

void ListesBranchesPourHGT(int *tab_sommets,long int * ARETE, int taille,struct DescTree *DT,int *tab_branches,int *nb_branches){
	
	(*nb_branches) = 0;
	for(int i=1;i<=2*taille-3;i++){
		DT[ARETE[2*i-1]].Tableau[0] = DT[ARETE[2*i-1]].nbSommet;
		DT[ARETE[2*i-2]].Tableau[0] = DT[ARETE[2*i-2]].nbSommet;
		
		if( ((Est_un_sous_ensemble_exact(DT[ARETE[2*i-1]].Tableau,tab_sommets)==1) || (Est_un_sous_ensemble_exact(DT[ARETE[2*i-2]].Tableau,tab_sommets)==1)) ||
		    ((Est_un_sous_ensemble_exact(DT[ARETE[2*i-2]].Tableau,tab_sommets)==1) || (Est_un_sous_ensemble_exact(DT[ARETE[2*i-1]].Tableau,tab_sommets)==1)) ){
			(*nb_branches) ++;
			
			tab_branches[(*nb_branches)] = i;
					
		}	
	}
	
}


void SAVEASNewick(double *LONGUEUR, long int *ARETE, char ** noms, int nn,const char* t){
	int n,root,a;
	int Ns;
	int i, j, sd, sf, *Suc, *Fre, *Tree, *degre, *Mark;
	double *Long;
	int *boot; 
	char *string = (char*)malloc(10000);	
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
	
	for (i=1;i<=2*n-3;i++){ 
		
		if (i<=n) degre[i]=1;
		else degre[i]=3;
	} 
	
	degre[2*n-2]=3;
	root=Ns+1;
	
	int cpt=0;
	
	for (;;){
		cpt++;
		if(cpt > 1000000) exit(1);
		a=0; a++;
		for (j=1;j<=2*n-2;j++)
			Mark[j]=0;
		
		for (i=1;i<=2*n-3;i++){ 	  									
			if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0)){
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-2]] = (int) bootStrap[i-1];
				
			}else if ((degre[ARETE[2*i-1]]==1)&&(degre[ARETE[2*i-2]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0)){
				Tree[ARETE[2*i-1]]=ARETE[2*i-2]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-1]]=LONGUEUR[i-1];
				
				if(bootStrap != NULL) boot[ARETE[2*i-1]] = (int) bootStrap[i-1];
				
			}else if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]==1)&&(Mark[ARETE[2*i-2]]==0)&&(Mark[ARETE[2*i-1]]==0)){
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; root=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; a=-1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				
				if(bootStrap != NULL) boot[ARETE[2*i-2]] = (int) bootStrap[i-1];
			}
			if (a==-1) break;
		}
		if (a==-1) break;
	}
	
	
	/*  On decale et on complete la structure d'arbre avec Successeurs et Freres  */
	for (i=Ns+1;i>0;i--){ 	
		Fre[i]=0; Suc[i]=0;
	}	
	
	Tree[root]=0;/*Tree[Ns+1]=0;*/
	
	for (i=1;i<=Ns+1/*Ns*/;i++){	
		if (i!=root){
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
	for (;;){	
		if(cpt > 1000000) exit(1);
		
		if (Suc[i]>0){	
			sprintf(string,"%s(",string);
			Suc[i]=-Suc[i]; i=-Suc[i]; }
		else if (Fre[i]!=0){	
			if (Suc[i]==0) sprintf(string,"%s%s:%.4f,",string,noms[i],Long[i]);
			else {
				if(bootStrap != NULL)
					sprintf(string,"%s%d:%.4f,",string,boot[i],Long[i]);
				else
					sprintf(string,"%s:%.4f,",string,Long[i]);
			}
			i=Fre[i]; 
		}
		else if (Tree[i]!=0){	
			if (Suc[i]==0) sprintf(string,"%s%s:%.4f)",string,noms[i],Long[i]);
		else {
			if(bootStrap != NULL)
				sprintf(string,"%s%d:%.4f)",string,boot[i],Long[i]);
			else
				sprintf(string,"%s:%.4f)",string,Long[i]);
			}
			i=Tree[i]; 
		}
		else break;
	}	
	strcat(string,";");
	
	FILE *pt_t = fopen(t,"w+");
	fprintf(pt_t,"%s",string);
	fclose(pt_t);
	
	free(Suc); free(Fre); free(Tree); free(Long); free(degre); free(Mark);	free(string);
}


void deleteBipartition(DescTree *DT,InputTree aTree){
	int i,j,k;
	for(i=1;i<=2*aTree.size-2-aTree.kt;i++){
		if(i!=aTree.size){
			for(j=0;j<=DT[i].nbSommet+1;j++)
				free(DT[i].Matrice[j]);
			free(DT[i].Matrice);
		}
		free(DT[i].Tableau);
	}
	free(DT);
}

void PrintHeader(FILE *out,Parameters param){
	char criterion[50];
	if(strcmp(param.criterion,"rf")==0)
		strcpy(criterion,"Robinson and Foulds distance");
	if(strcmp(param.criterion,"ls")==0)
		strcpy(criterion,"least-squares");
	if(strcmp(param.criterion,"bd")==0)
		strcpy(criterion,"bipartition dissimilarity");
	if(strcmp(param.criterion,"qd")==0)
		strcpy(criterion,"quartet distance");

	fprintf(out,"%s",description);
	fprintf(out,"\n");
	fprintf(out,"\nSubtree constraint   :%s",param.subtree);
	fprintf(out,"\nCriterion            :%s",criterion);
	fprintf(out,"\nBootstrap            :%s",param.bootstrap);
	fprintf(out,"\n");
}

//=================================================================
//== 
//=================================================================
void initInputTree(struct InputTree *aTree){
	aTree->Adjacence = NULL;
	aTree->ARETE = NULL;
	aTree->LONGUEUR = NULL;
	aTree->ADD = NULL;
	aTree->Root = -1;
	aTree->size = -1;
	aTree->SpeciesName = NULL;
	aTree->Input = NULL;
	aTree->kt = 0;
}

//============================================
//=
//============================================
void printRoot(char *fichier, int R1, int R2){
	
	FILE *out;
	
	if((out=fopen(fichier,"w+"))==NULL){
		printf("\nCan't open root file (%s)",fichier);
		exit(-1);
	}
	fprintf(out,"%d %d",R1,R2);
	fclose(out);
}

void printRootByLeaves(char *fichier,int choix,struct InputTree *aTree){

	FILE *out;
	int i;
	if((out=fopen(fichier,"w+"))==NULL){
		printf("\nCan't open root file (%s)",fichier);
		exit(-1);
	}
	
	for(i=1;i<=aTree->size;i++)
		if(aTree->ADD[i][aTree->ARETE[2*choix-1]] < aTree->ADD[i][aTree->ARETE[2*choix-2]])
			fprintf(out,"%s ",aTree->SpeciesName[i]);
	
	fprintf(out,"\n<>\n");
	
	for(i=1;i<=aTree->size;i++)
		if(aTree->ADD[i][aTree->ARETE[2*choix-1]] > aTree->ADD[i][aTree->ARETE[2*choix-2]])
			fprintf(out,"%s ",aTree->SpeciesName[i]);
	
	fclose(out);
	
}
//=====================================================
//
//=====================================================
void allocMemmory(struct InputTree *aTree, int n){
	int i;

	if(aTree->ARETE == NULL){
		aTree->degre = (int*)malloc(2*n*sizeof(int));
		aTree->ADD = (double**)malloc(2*n*sizeof(double*));
		aTree->Adjacence = (double**)malloc(2*n*sizeof(double*));
		aTree->Input = (double**)malloc(2*n*sizeof(double*));
		aTree->W = (double**)malloc((n+1)*sizeof(double*));
		
		for(i=0;i<2*n;i++){
			aTree->ADD[i] = (double*)malloc(2*n*sizeof(double));
			aTree->Adjacence[i] = (double*)malloc(2*n*sizeof(double));
			aTree->Input[i] = (double*)malloc(2*n*sizeof(double));
			/* if(i<=n)
				aTree->W[i] = (double*)malloc(2*n*sizeof(double)); */
		}

		for(i=0;i<=n;i++)
			aTree->W[i] = (double*)malloc(2*n*sizeof(double));
		
		aTree->ARETE    =(long int*)malloc(4*(2*(n))*sizeof(long int));
		aTree->LONGUEUR	=(double*)malloc((4*(n))*sizeof(double));
		aTree->SpeciesName = (char**)malloc(2*n*sizeof(char*));
		
		for(i=0;i<=n;i++)
			aTree->SpeciesName[i] = (char*)malloc(50);
	}

}

void freeInputTree(struct InputTree *aTree,int n){
	int i,j;
	
	for(i=0;i<2*n;i++){
		free(aTree->ADD[i]);
		free(aTree->Adjacence[i]);
		free(aTree->Input[i]);
		if(i<=n){
			free(aTree->W[i]);
		}
	}
	
	for(i=0;i<=n;i++){
		free(aTree->SpeciesName[i]);
	}
	
	free(aTree->degre);
	free(aTree->ADD);
	free(aTree->Adjacence);
	free(aTree->Input);
	free(aTree->W);
	free(aTree->ARETE);
	free(aTree->LONGUEUR);
	free(aTree->SpeciesName);

}

void freeReducedTree(struct InputTree *aTree,int n){
	int i,j;
	int inc = 10;
	
	for(i=0;i<2*n;i++){
		free(aTree->ADD[i]);
		free(aTree->Input[i]);
		if(i<=n)
			free(aTree->W[i]);
	}
	
	for(i=0;i<2*(n+inc);i++){
		free(aTree->Adjacence[i]);
	}
	
	for(i=0;i<=n;i++){
		free(aTree->SpeciesName[i]);
	}
	
	free(aTree->degre);
	free(aTree->ADD);
	free(aTree->Adjacence);
	free(aTree->Input);
	free(aTree->W);
	free(aTree->ARETE);
	free(aTree->LONGUEUR);
	free(aTree->SpeciesName);

}


void freeTree(struct InputTree *aTree,int n){
	int i,j;
	int inc = 10;
	
	for(i=0;i<2*n;i++){
		free(aTree->ADD[i]);
		free(aTree->Input[i]);
		if(i<=n)
			free(aTree->W[i]);
	}
	
	for(i=0;i<=n;i++){
		free(aTree->SpeciesName[i]);
	}
	
	free(aTree->ADD);
	free(aTree->Input);
	free(aTree->W);
	free(aTree->SpeciesName);

}
//================================================================================
//==
//================================================================================

int TestSubTreeConstraint(struct InputTree aTree,int source, int dest,struct DescTree *DTSpecies, struct DescTree *DTGene){

	int basSource,basDest;
	int etape1=0,etape2=0;
	int LGTpossible = 0;
	int pv,i,y,z,tk,tl;
	int n = aTree.size;
	int mk1,mk2,nbSom;
	struct DescTree DTemp;
	double ** DISTemp =(double **)malloc((2*n+1)*sizeof(double*));

	int * PLACEk1=(int *) malloc((2*n-2)*sizeof(int));
	int * PLACEk2=(int *) malloc((2*n-2)*sizeof(int));
	
	int ** Bk1=(int **) malloc((2*n-2)*sizeof(int*));
	int ** Bk2=(int **) malloc((2*n-2)*sizeof(int*));

	for (i=0;i<2*n-2;i++)
	{
		Bk1[i]=(int *) malloc((n)*sizeof(int));
		Bk2[i]=(int *) malloc((n)*sizeof(int));     
		DISTemp[i]=(double *)malloc((2*n+1)*sizeof(double));  
	}
	
	//== find the node inf for the source branch
	if(aTree.ADD[aTree.ARETE[2*source-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*source-2]][aTree.Root])
		basSource = aTree.ARETE[2*source-2];
	else
		basSource = aTree.ARETE[2*source-1];

	//== find the node inf for the dest branch
	if(aTree.ADD[aTree.ARETE[2*dest-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*dest-2]][aTree.Root])
		basDest = aTree.ARETE[2*dest-2];
	else
		basDest = aTree.ARETE[2*dest-1];


	if(basSource<n)
		etape1=1;
	else{
		for(pv=n+1;pv<2*n-2;pv++){
			if(vecteursEgaux(DTSpecies[basSource],DTGene[pv]) == 1){	
				/*if(DTSpecies[basSource].nbSommet == 3)
					etape1= 0;
				else*/{
					mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet+1);
					mk2=Bipartition_Table(DTSpecies[basSource].Matrice,Bk2,PLACEk2,DTSpecies[basSource].nbSommet+1);
					if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTSpecies[basSource].nbSommet+1) == 0)
						etape1 = 1;
				}
				pv = (int)INFINI;
			}
		}
	}

	if(etape1){/*/== comparer les sous arbres au niveau destination*/
		if(basDest<n)
			etape2 = 1;
		else{
			for(pv=n+1;pv<2*n-2;pv++){
				if(vecteursEgaux(DTSpecies[basDest],DTGene[pv]) == 1){
					/*if(DTSpecies[basDest].nbSommet == 3)
						etape2 = 0;
					else*/{
						mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet+1);
						mk2=Bipartition_Table(DTSpecies[basDest].Matrice,Bk2,PLACEk2,DTSpecies[basDest].nbSommet+1);
						if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTGene[pv].nbSommet+1) == 0)
							etape2 = 1;
					}
					pv = (int)INFINI;
				}				
			}	
		}
	}

	if(etape1==1 && etape2==1){
		
		nbSom = DTemp.nbSommet = DTSpecies[basDest].nbSommet + DTSpecies[basSource].nbSommet;

		DTemp.Tableau = (int*)malloc((DTemp.nbSommet+1)*sizeof(int));

		for(y=1;y<=DTSpecies[basDest].nbSommet;y++)
			DTemp.Tableau[y] = DTSpecies[basDest].Tableau[y];
		y--;
		for(z=1;z<=DTSpecies[basSource].nbSommet;z++)
			DTemp.Tableau[z+y] = DTSpecies[basSource].Tableau[z];

		TrierTableau(DTemp.Tableau,DTemp.nbSommet);

		for(pv=n+1;pv<=2*n-2;pv++){
			if(vecteursEgaux(DTemp,DTGene[pv])){

				for(tk=1;tk<=DTGene[pv].nbSommet;tk++){
					for(tl=1;tl<=DTGene[pv].nbSommet;tl++)
						DISTemp[tk][tl] = aTree.ADD[DTemp.Tableau[tk]][DTemp.Tableau[tl]];
				}
				
				mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet);
				mk2=Bipartition_Table(DISTemp,Bk2,PLACEk2,DTemp.nbSommet);
				if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTGene[pv].nbSommet) == 0)
					LGTpossible = 1;
			}
		}
		free(DTemp.Tableau);	
	}

	free(PLACEk1);
	free(PLACEk2);

	for (i=0;i<2*n-2;i++)
	{
		free(Bk1[i]);
		free(Bk2[i]);
		free(DISTemp[i]);
	}
	free(DISTemp);	
	free(Bk1);
	free(Bk2);

	return LGTpossible;
		
}

//================================================================================
//==
//================================================================================

int TestSubTreeLeafs(struct InputTree aTree,int source, int dest,struct DescTree *DTSpecies, struct DescTree *DTGene){

	int basSource,hautSource,basDest,noeudRacine;
	int etape1=0,etape2=0;
	int LGTpossible = 0;
	int pv,i,y,z,tk,tl;
	int n = aTree.size;
	int mk1,mk2,nbSom;
	struct DescTree DTemp;
	double ** DISTemp =(double **)malloc((2*n+1)*sizeof(double*));

	int * PLACEk1=(int *) malloc((2*n-2)*sizeof(int));
	int * PLACEk2=(int *) malloc((2*n-2)*sizeof(int));
	
	int ** Bk1=(int **) malloc((2*n-2)*sizeof(int*));
	int ** Bk2=(int **) malloc((2*n-2)*sizeof(int*));

	for (i=0;i<2*n-2;i++)
	{
		Bk1[i]=(int *) malloc((n)*sizeof(int));
		Bk2[i]=(int *) malloc((n)*sizeof(int));     
		DISTemp[i]=(double *)malloc((2*n+1)*sizeof(double));  
	}
	
	for (i=1;i<2*n-3-aTree.kt;i++){
		if(aTree.ARETE[2*i-1] == aTree.Root)
			noeudRacine = aTree.ARETE[2*i-2];
		if(aTree.ARETE[2*i-2] == aTree.Root)
			noeudRacine = aTree.ARETE[2*i-1];
	}
			
	//== find the node inf for the source branch
	if(aTree.ADD[aTree.ARETE[2*source-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*source-2]][aTree.Root]){
		basSource = aTree.ARETE[2*source-2];
		hautSource = aTree.ARETE[2*source-1];
	}
	else{
		basSource = aTree.ARETE[2*source-1];
		hautSource = aTree.ARETE[2*source-2];
	}
	//== find the node inf for the dest branch
	if(aTree.ADD[aTree.ARETE[2*dest-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*dest-2]][aTree.Root])
		basDest = aTree.ARETE[2*dest-2];
	else
		basDest = aTree.ARETE[2*dest-1];
		
		
	if(basSource<n)
		etape1=1;
	else{
		for(pv=n+1;pv<2*n-2;pv++){
			if(vecteursEgaux(DTSpecies[basSource],DTGene[pv]) == 1){	
				/*if(DTSpecies[basSource].nbSommet == 3)
					etape1= 0;
				else*/{
					mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet+1);
					mk2=Bipartition_Table(DTSpecies[basSource].Matrice,Bk2,PLACEk2,DTSpecies[basSource].nbSommet+1);
					if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTSpecies[basSource].nbSommet+1) == 0)
						etape1 = 1;
				}
				pv = (int)INFINI;
			}
		}
	}

	if(etape1){/*/== comparer les sous arbres au niveau destination*/
		if(basDest<n)
			etape2 = 1;
		else{
			for(pv=n+1;pv<2*n-2;pv++){
				if(vecteursEgaux(DTSpecies[basDest],DTGene[pv]) == 1){
					/*if(DTSpecies[basDest].nbSommet == 3)
						etape2 = 0;
					else*/{
						mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet+1);
						mk2=Bipartition_Table(DTSpecies[basDest].Matrice,Bk2,PLACEk2,DTSpecies[basDest].nbSommet+1);
						if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTGene[pv].nbSommet+1) == 0)
							etape2 = 1;
					}
					pv = (int)INFINI;
				}				
			}	
		}
	}

	if(etape1==1 && etape2==1){
		
		nbSom = DTemp.nbSommet = DTSpecies[basDest].nbSommet + DTSpecies[hautSource].nbSommet;
				
		DTemp.Tableau = (int*)malloc((DTemp.nbSommet+1)*sizeof(int));
		
		for(y=1;y<=DTSpecies[basDest].nbSommet;y++){
			DTemp.Tableau[y] = DTSpecies[basDest].Tableau[y];
		}
		
		y--;
		
		for(z=1;z<=DTSpecies[hautSource].nbSommet;z++){
			DTemp.Tableau[z+y] = DTSpecies[hautSource].Tableau[z];
		}
	
		
		TrierTableau(DTemp.Tableau,DTemp.nbSommet);

		for(pv=n+1;pv<=2*n-2;pv++){
			if(vecteursEgaux(DTemp,DTGene[pv])){
				LGTpossible = 1;
			}
		}
		
		for(pv=n+1;pv<=2*n-2;pv++){
			if(vecteursEgaux(DTSpecies[hautSource],DTGene[pv])){
				LGTpossible = 1;
			}
		}
		
		
		free(DTemp.Tableau);	
	}

	free(PLACEk1);
	free(PLACEk2);

	for (i=0;i<2*n-2;i++)
	{
		free(Bk1[i]);
		free(Bk2[i]);
		free(DISTemp[i]);
	}
	free(DISTemp);	
	free(Bk1);
	free(Bk2);
	
	if(hautSource == noeudRacine)
		LGTpossible = 1;
	
	return LGTpossible;
		
}
//=================================================================
//==
//=================================================================
void loadCriteria(struct CRITERIA aCrit, struct HGT *aHGT){

	aHGT->crit.LS = aCrit.LS;
	aHGT->crit.RF = aCrit.RF;
	aHGT->crit.BD = aCrit.BD;
	aHGT->crit.QD = aCrit.QD;
	aHGT->valide = 1;
}

//=================================================================
//==
//=================================================================
void computeCriteria(double ** Matrix1, double ** Matrix2, int size,struct CRITERIA *aCrit,double *L1, long int *A1,double *L2, long int *A2){

	int mI,m,i,j,RF,QD;
	double LS,BD=0;

	//= robinson and foulds
	m=Bipartition_Table(Matrix1,aCrit->B,aCrit->PLACE,size);
	mI=Bipartition_Table(Matrix2,aCrit->BI,aCrit->PLACEI,size);
	RF = Table_Comparaison(aCrit->B,aCrit->BI,aCrit->PLACE,aCrit->PLACEI,m,mI,size); 

	//= least-squares
	LS = 0.0;
	for (i=1;i<=size-1;i++)
	{
		for (j=i+1;j<=size;j++){
			LS=LS + (Matrix1[i][j]-Matrix2[i][j])*(Matrix1[i][j]-Matrix2[i][j]);
			if(LS > INFINI){
				
			}
		}
	}
	//= Bipartition Distance
	BD = BipartitionDistance(aCrit->B,aCrit->BI,size);
	
	//= Quartet Distance
	//= le calcul de QD necessite la creation de 2 fichiers d'arbres t1 et t2
	//= sinon la distance sera de -1 (distance non calculee)
	
	aCrit->LS = LS;
	aCrit->BD = BD;
	aCrit->RF = RF;
	aCrit->QD = QD;

}

double MIN2 (double A,double B){
	if ( A > B) return B;
	return A;
}
//=================================================================
//==
//=================================================================
void ComputeNewDistances(double **DIST,double **ADD,int size,int kt,int newNode,double tailleBranche,double tailleTransfert,int Ssup,int Sinf,int Sdest){

	
	
	int i,j;
	double lt = DIST[Sdest][size];


	for(i=1;i<=2*size-2-kt;i++){
		
		for(j=i+1;j<=2*size-2-kt;j++){
			
			if(i!=newNode && j!=newNode){
				if( (fabs(DIST[i][size]-DIST[i][Sdest]-lt)< 3*epsilon) && (fabs(DIST[j][size]-DIST[j][Sdest]-lt)>3*epsilon) ){
				  if((DIST[i][Sdest]+tailleTransfert+tailleBranche+DIST[Ssup][j]) < (DIST[i][Sdest]+tailleTransfert+tailleBranche+DIST[Sinf][j])){
                    ADD[i][j]  = ADD[j][i] = DIST[i][Sdest]+tailleTransfert+tailleBranche+DIST[Ssup][j];
                  }else{
                    ADD[i][j]  = ADD[j][i] = DIST[i][Sdest]+tailleTransfert+tailleBranche+DIST[Sinf][j];
                  }
				}else if( (fabs(DIST[i][size]-DIST[i][Sdest]-lt)>3*epsilon) && (fabs(DIST[j][size]-DIST[j][Sdest]-lt)<3*epsilon) ){
					if((DIST[j][Sdest]+tailleTransfert+tailleBranche+DIST[Ssup][i]) < (DIST[j][Sdest]+tailleTransfert+tailleBranche+DIST[Sinf][i])){
					   ADD[i][j] = ADD[j][i] = DIST[j][Sdest]+tailleTransfert+tailleBranche+DIST[Ssup][i]; 
					}else{
                        ADD[i][j] = ADD[j][i] = DIST[j][Sdest]+tailleTransfert+tailleBranche+DIST[Sinf][i];
                    }
                }
            }

		}

	}

}

//=================================================================
//==
//=================================================================
void applyHGT(double**ref,struct InputTree * aTree,int i,int j){
	
	int nodeToDel,otherNode,neighbor1=0,neighbor2=0,p,q,branch1=0,branch2=0,newNode;
	int nouveauNoeud;
	double tailleBrancheSource,tailleTransfert;
	int Ssup,Sinf,Sdestination;
	int kt = aTree->kt;
	int iternumber=100;

	if(i==0 || j==0) { printf("\ntransfert impossible"); exit(-1);}
	double ** DIST = (double**)malloc(2*aTree->size*(sizeof(double*)));
	for(p=0;p<2*aTree->size;p++)
		DIST[p] = (double *)malloc(2*aTree->size*(sizeof(double)));

  //= copie du contenue de ADD dans DIST
	for(p=1;p<=2*aTree->size-2-kt;p++)	
		for(q=1;q<=2*aTree->size-2-kt;q++)
			DIST[p][q] = aTree->ADD[p][q];
	
	//= recherche des sommets superieurs et inferieurs sources
	tailleBrancheSource = aTree->LONGUEUR[i];
	if ( aTree->ADD[aTree->ARETE[2*i-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*i-2]][aTree->Root]){
			Ssup = aTree->ARETE[2*i-1];
			Sinf = aTree->ARETE[2*i-2];
	}
	else{
			Ssup = aTree->ARETE[2*i-2];
			Sinf = aTree->ARETE[2*i-1];
	}
	
	//= recherche des sommets superieurs et inferieurs destinations
	if(aTree->ADD[aTree->ARETE[2*j-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*j-2]][aTree->Root])
	{nodeToDel = aTree->ARETE[2*j-1]; otherNode = Sdestination = aTree->ARETE[2*j-2];}
	else
	{nodeToDel = aTree->ARETE[2*j-2]; otherNode = Sdestination = aTree->ARETE[2*j-1];}

	//= cas d'un noeud binaire, on va deplacer le noeud superieur (nodeToDel)
	if(aTree->degre[nodeToDel] == 3){
		//== 1. recherche des voisins du noeud a supprimer
		for(p=1;p<=2*(aTree->size)-3-aTree->kt;p++){
			if(aTree->ARETE[2*p-1] == nodeToDel && aTree->ARETE[2*p-2] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-2]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-2]; branch2=p;}
			} 
			if(aTree->ARETE[2*p-2] == nodeToDel && aTree->ARETE[2*p-1] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-1]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-1]; branch2=p;}
			} 
		}
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) ){
					  if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) <(DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){ 
						  aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
					  }else{
                          aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
					  }
					}
				}
		}
		
		//== 3. branch1 is used to connect the two neighbors
		aTree->ARETE[2*branch1-1] = neighbor1;
		aTree->ARETE[2*branch1-2] = neighbor2;
		aTree->LONGUEUR[branch1-1]= aTree->LONGUEUR[branch1-1] + aTree->LONGUEUR[branch2-1];

		//== 4. branch2 is used to connect one of the source node (i)
		aTree->ARETE[2*branch2-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*branch2-2] = nodeToDel;
		aTree->LONGUEUR[branch2-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;

		//== 5. branch i is used to connect the other source node (i)
		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = nodeToDel;
		aTree->LONGUEUR[i-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		nouveauNoeud = nodeToDel;

		for(p=1;p<=2*aTree->size-2-aTree->kt;p++){
			
			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
				aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = 1.0 + DIST[Sdestination][p];
			}else{
			  if(DIST[p][Ssup] < DIST[p][Sinf]){
					aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Ssup] + aTree->LONGUEUR[i-1];
			  }else{
					aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Sinf] + aTree->LONGUEUR[i-1];
			  }
		    }
		}
		
		aTree->ADD[nodeToDel][nodeToDel] = 0.0;
		tailleBrancheSource = aTree->LONGUEUR[i-1];
		tailleTransfert = 1.0; 
		aTree->LONGUEUR[j-1]=1.0;
		nouveauNoeud=nodeToDel;
	}
	//= noeud de degre > a 3
	else{
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					   if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) < (DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){
						    aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						 }
						 else{
                aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
             }
					}
				}
		}
		
		//== 2. Ajouter un nouveau noeud sur la branche i (ajout d'une nouvelle arete)
		newNode = nouveauNoeud = 2*aTree->size-2-aTree->kt+1;

		double lt = (aTree->LONGUEUR[i-1]/2.0 >= 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		//== 4. connecter le nouveau noeud
		aTree->ARETE[2*j-1]=newNode;
		aTree->ARETE[2*j-2]=otherNode;
		aTree->LONGUEUR[j-1]= 1.0;

		//== 5. connecter la nouvelle arete
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-2] = newNode;
		aTree->LONGUEUR[(2*aTree->size-3-aTree->kt+1)-1]= lt; //1.0;

		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = newNode;
		aTree->LONGUEUR[i-1]= lt; //1.0;

		aTree->degre[newNode] = 3;
		aTree->degre[nodeToDel]--;
		aTree->kt--;
		
		for(p=1;p<=newNode;p++){
			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
				aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = 1.0 + DIST[Sdestination][p];
			}
			else	{
				if((DIST[p][Ssup] + lt) < (DIST[p][Sinf] + lt)){
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Ssup] + lt;
        }
        else{
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Sinf] + lt;
        }
			}
		}	
	
		
		
		aTree->ADD[newNode][newNode] = 0.0;
		tailleBrancheSource = lt;
		tailleTransfert=1.0;
	}
	
	//== Mise a jour des distances
		
	int nbNoeud1=0;
	int nbNoeud2=0;
	int * tab_noeud1 = (int*)malloc((2*aTree->size+1)*sizeof(int));
	int * tab_noeud2 = (int*)malloc((2*aTree->size+1)*sizeof(int));
	
  //== liste des noeuds sous sDestination
  
	for(int i=1;i<2*aTree->size-2-aTree->kt;i++){
      if( fabs(DIST[i][aTree->size]-DIST[i][Sdestination] - DIST[Sdestination][aTree->size]) < 5*epsilon ){
        tab_noeud1[nbNoeud1++] = i;
      }
      else{
        tab_noeud2[nbNoeud2++] = i;
      }
  }  
    //== mise a jour des distances
  	for(int i=0;i<nbNoeud1;i++){
      for(int j=0;j<nbNoeud2;j++){
        aTree->ADD[tab_noeud1[i]][tab_noeud2[j]] =  aTree->ADD[tab_noeud2[j]][tab_noeud1[i]] = aTree->ADD[tab_noeud2[j]][nouveauNoeud] + tailleTransfert + aTree->ADD[Sdestination ][tab_noeud1[i]];
      }
    }
  	
	kt = aTree->kt;
	
	loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size, aTree->kt);
	
	for(p=0;p<2*aTree->size;p++)
		free(DIST[p]);
	free(DIST);
}

//=================================================================
//==
//=================================================================
void applyHGT2(double**ref,struct InputTree * aTree,int i,int j){
	
	int nodeToDel,otherNode,neighbor1=0,neighbor2=0,p,q,branch1=0,branch2=0,newNode;
	int nouveauNoeud;
	double tailleBrancheSource,tailleTransfert;
	int Ssup,Sinf,Sdestination;
	int kt = aTree->kt;
	int iternumber=100;

	if(i==0 || j==0) { printf("\ntransfert impossible"); exit(-1);}
	double ** DIST = (double**)malloc(2*aTree->size*(sizeof(double*)));
	for(p=0;p<2*aTree->size;p++)
		DIST[p] = (double *)malloc(2*aTree->size*(sizeof(double)));

  //= copie du contenue de ADD dans DIST
	for(p=1;p<=2*aTree->size-2-kt;p++)	
		for(q=1;q<=2*aTree->size-2-kt;q++)
			DIST[p][q] = aTree->ADD[p][q];
	
	//= recherche des sommets superieurs et inferieurs sources
	tailleBrancheSource = aTree->LONGUEUR[i];
	if ( aTree->ADD[aTree->ARETE[2*i-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*i-2]][aTree->Root]){
			Ssup = aTree->ARETE[2*i-1];
			Sinf = aTree->ARETE[2*i-2];
	}
	else{
			Ssup = aTree->ARETE[2*i-2];
			Sinf = aTree->ARETE[2*i-1];
	}
	
	//= recherche des sommets superieurs et inferieurs destinations
	if(aTree->ADD[aTree->ARETE[2*j-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*j-2]][aTree->Root])
	{nodeToDel = aTree->ARETE[2*j-1]; otherNode = Sdestination = aTree->ARETE[2*j-2];}
	else
	{nodeToDel = aTree->ARETE[2*j-2]; otherNode = Sdestination = aTree->ARETE[2*j-1];}

	//= cas d'un noeud binaire, on va deplacer le noeud superieur (nodeToDel)
	if(aTree->degre[nodeToDel] == 3){
		//== 1. recherche des voisins du noeud a supprimer
		for(p=1;p<=2*(aTree->size)-3-aTree->kt;p++){
			if(aTree->ARETE[2*p-1] == nodeToDel && aTree->ARETE[2*p-2] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-2]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-2]; branch2=p;}
			} 
			if(aTree->ARETE[2*p-2] == nodeToDel && aTree->ARETE[2*p-1] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-1]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-1]; branch2=p;}
			} 
		}
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					  if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) <(DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){ 
						  aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						}
						else{
               aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
            }
					}
				}
		}
		
		//== 3. branch1 is used to connect the two neighbors
		aTree->ARETE[2*branch1-1] = neighbor1;
		aTree->ARETE[2*branch1-2] = neighbor2;
		aTree->LONGUEUR[branch1-1]= aTree->LONGUEUR[branch1-1] + aTree->LONGUEUR[branch2-1];

		//== 4. branch2 is used to connect one of the source node (i)
		aTree->ARETE[2*branch2-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*branch2-2] = nodeToDel;
		aTree->LONGUEUR[branch2-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;

		//== 5. branch i is used to connect the other source node (i)
		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = nodeToDel;
		aTree->LONGUEUR[i-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		nouveauNoeud = nodeToDel;

		for(p=1;p<=2*aTree->size-2-aTree->kt;p++){
			
			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
			 //printf("\n==%d,%d,%d",p,Sdestination,nodeToDel);
      	aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = 1.0 + DIST[Sdestination][p];
			}
			else	{
			  if(DIST[p][Ssup] < DIST[p][Sinf]){
				  aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Ssup] + aTree->LONGUEUR[i-1];
				}
				else{
          aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Sinf] + aTree->LONGUEUR[i-1];
        }
			}
		}
		
		aTree->ADD[nodeToDel][nodeToDel] = 0.0;
		tailleBrancheSource = aTree->LONGUEUR[i-1];
		tailleTransfert = 1.0; 
		aTree->LONGUEUR[j-1]=1.0;
		nouveauNoeud=nodeToDel;
	}
	//= noeud de degre > a 3
	else{
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					   if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) < (DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){
						    aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						 }
						 else{
                aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
             }
					}
				}
		}
		
		//== 2. Ajouter un nouveau noeud sur la branche i (ajout d'une nouvelle arete)
		newNode = nouveauNoeud = 2*aTree->size-2-aTree->kt+1;

		double lt = (aTree->LONGUEUR[i-1]/2.0 >= 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		//== 4. connecter le nouveau noeud
		aTree->ARETE[2*j-1]=newNode;
		aTree->ARETE[2*j-2]=otherNode;
		aTree->LONGUEUR[j-1]= 1.0;

		//== 5. connecter la nouvelle arete
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-2] = newNode;
		aTree->LONGUEUR[(2*aTree->size-3-aTree->kt+1)-1]= lt; //1.0;

		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = newNode;
		aTree->LONGUEUR[i-1]= lt; //1.0;

		aTree->degre[newNode] = 3;
		aTree->degre[nodeToDel]--;
		aTree->kt--;
		
		for(p=1;p<=newNode;p++){
			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
				aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = 1.0 + DIST[Sdestination][p];
			}
			else	{
				if((DIST[p][Ssup] + lt) < (DIST[p][Sinf] + lt)){
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Ssup] + lt;
        }
        else{
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Sinf] + lt;
        }
			}
		}	
		
		aTree->ADD[newNode][newNode] = 0.0;
		tailleBrancheSource = lt;
		tailleTransfert=1.0;
	}
	
	int nbNoeud1=0;
	int nbNoeud2=0;
	int * tab_noeud1 = (int*)malloc((2*aTree->size+1)*sizeof(int));
	int * tab_noeud2 = (int*)malloc((2*aTree->size+1)*sizeof(int));
	
  //== liste des noeuds sous sDestination
  
	for(int i=1;i<2*aTree->size-2-aTree->kt;i++){
      if( fabs(DIST[i][aTree->size]-DIST[i][Sdestination] - DIST[Sdestination][aTree->size]) < 5*epsilon ){
        tab_noeud1[nbNoeud1++] = i;
      }
      else{
        tab_noeud2[nbNoeud2++] = i;
      }
  }  
    //== mise a jour des distances
  	for(int i=0;i<nbNoeud1;i++){
      for(int j=0;j<nbNoeud2;j++){
        aTree->ADD[tab_noeud1[i]][tab_noeud2[j]] =  aTree->ADD[tab_noeud2[j]][tab_noeud1[i]] = aTree->ADD[tab_noeud2[j]][nouveauNoeud] + tailleTransfert + aTree->ADD[Sdestination ][tab_noeud1[i]];
      }
    }
  	
	
	kt = aTree->kt;
	
  approx_arb(ref,aTree->ADD,aTree->ADD,aTree->W,&iternumber,aTree->ARETE,aTree->LONGUEUR,1,&kt,0,aTree->size);
  loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size, aTree->kt);
  Floyd(aTree->Adjacence,aTree->ADD,aTree->size,aTree->kt);

	for(p=0;p<2*aTree->size;p++)
		free(DIST[p]);
	free(DIST);
}
//=================================================================
//==
//=================================================================
void findListSpecies(struct HGT *bestHGT, struct DescTree *DTSpecies,struct InputTree aTree){

	int source,dest,i;
	
	if(aTree.ADD[aTree.ARETE[2*(bestHGT->source)-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*(bestHGT->source)-2]][aTree.Root])
		source = aTree.ARETE[2*(bestHGT->source)-2];
	else
		source = aTree.ARETE[2*(bestHGT->source)-1];

	//== find the node inf for the dest branch
	if(aTree.ADD[aTree.ARETE[2*(bestHGT->destination)-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*(bestHGT->destination)-2]][aTree.Root])
		dest = aTree.ARETE[2*(bestHGT->destination)-2];
	else
		dest = aTree.ARETE[2*(bestHGT->destination)-1];
	
	if(bestHGT->listSource != NULL)
		free(bestHGT->listSource);
	bestHGT->listSource = (int *)malloc((DTSpecies[source].nbSommet+1)*sizeof(int));
	bestHGT->listSource[0] = DTSpecies[source].nbSommet;
	for(i=1;i<=DTSpecies[source].nbSommet;i++) bestHGT->listSource[i] = DTSpecies[source].Tableau[i];

	if(bestHGT->listDestination != NULL)
		free(bestHGT->listDestination);
	bestHGT->listDestination = (int *)malloc((DTSpecies[dest].nbSommet+1)*sizeof(int));
	bestHGT->listDestination[0] = DTSpecies[dest].nbSommet;
	for(i=1;i<=DTSpecies[dest].nbSommet;i++) bestHGT->listDestination[i] = DTSpecies[dest].Tableau[i];
}

//=================================================================
//==
//=================================================================
int TestCriterionAndUpdate(int *first,const char *CRITERION, struct CRITERIA aCrit, struct HGT * aHGT,int i,int j,int flag,int flag_bootstrap){

	int best = 0;

	if (strcmp(CRITERION,"rf")==0){
		if (aHGT->crit.RF > aCrit.RF) {
			best = 1;			
			(*first)=0;
		}
		if (aHGT->crit.RF == aCrit.RF && aHGT->crit.LS > aCrit.LS && (*first) == 0)
			best = 1; 
	}
	else if (strcmp(CRITERION,"ls")==0){
		if (aHGT->crit.LS > aCrit.LS) best = 1;
	}
	else if (strcmp(CRITERION,"bd")==0){
		if (aHGT->crit.BD > aCrit.BD) best = 1;
		if (flag == 2){
			printf("LALA");
			best=1;
		}
		
	}
	else if (strcmp(CRITERION,"qd")==0){
		if (aHGT->crit.QD > aCrit.QD) best = 1;
	}
	else{
		printf("\nunknown criterion %s, try again...\n",CRITERION);
		exit(-1);
	}

	if (flag_bootstrap == 1){
		
		if(rand_bootstrap == 0){
			best = 0;
		}
		else{
			best = 1;
		}
	}
	
	if(best == 1){
		aHGT->destination = j;
		aHGT->source = i;
		aHGT->crit.BD = aCrit.BD;
		aHGT->crit.LS = aCrit.LS;
		aHGT->crit.RF = aCrit.RF;
		aHGT->crit.QD = aCrit.QD;
	}
	
	return best;
}

//=================================================================
//==
//=================================================================
void UpdateCriterion(int *first,const char *CRITERION, struct CRITERIA aCrit, struct HGT * aHGT,int i,int j,int flag_bootstrap){

		aHGT->destination = j;
		aHGT->source = i;
		aHGT->crit.BD = aCrit.BD;
		aHGT->crit.LS = aCrit.LS;
		aHGT->crit.RF = aCrit.RF;
		aHGT->crit.QD = aCrit.QD;
}

//=================================================================
//==
//=================================================================
int isAValidHGT(struct InputTree SpeciesTree,int i, int j){
	
	int d1,d2,d3,d4,root;

	//== if is the root branch
	if(SpeciesTree.ARETE[2*i-1] == SpeciesTree.Root || SpeciesTree.ARETE[2*i-2] == SpeciesTree.Root || 
	   SpeciesTree.ARETE[2*j-1] == SpeciesTree.Root || SpeciesTree.ARETE[2*j-2] == SpeciesTree.Root )
		return 0;

	//== if there is a direct common ancestor
	if(((SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-1])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-1]]==3)) || 
	   ((SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-2])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-1]]==3)) ||
       ((SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-1])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-2]]==3)) || 
	   ((SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-2]]==3)) )
		return 0;

	//== if j is a branch connected to the root branch
	if(SpeciesTree.ARETE[2*j-1] == 2*(SpeciesTree.size)-2 || SpeciesTree.ARETE[2*j-2] == 2*(SpeciesTree.size)-2)
		return 0;

	//== if i and j are on the same lineage

	root = SpeciesTree.Root;

	if(SpeciesTree.ADD[root][SpeciesTree.ARETE[2*i-1]] < SpeciesTree.ADD[root][SpeciesTree.ARETE[2*i-2]]){
		d1 = SpeciesTree.ARETE[2*i-1];
		d2 = SpeciesTree.ARETE[2*i-2];
	}
	else{
		d1 = SpeciesTree.ARETE[2*i-2];
		d2 = SpeciesTree.ARETE[2*i-1];
	}

	if(SpeciesTree.ADD[root][SpeciesTree.ARETE[2*j-1]] < SpeciesTree.ADD[root][SpeciesTree.ARETE[2*j-2]]){
		d3 = SpeciesTree.ARETE[2*j-1];
		d4 = SpeciesTree.ARETE[2*j-2];
	}
	else{
		d3 = SpeciesTree.ARETE[2*j-2];
		d4 = SpeciesTree.ARETE[2*j-1];
	}
	
	if(fabs(SpeciesTree.ADD[root][d4]-SpeciesTree.ADD[root][d1]-SpeciesTree.ADD[d1][d2]-SpeciesTree.ADD[d2][d3]-SpeciesTree.ADD[d3][d4]) < epsilon)
		return 0;
	if(fabs(SpeciesTree.ADD[root][d2]-SpeciesTree.ADD[root][d3]-SpeciesTree.ADD[d3][d4]-SpeciesTree.ADD[d4][d1]-SpeciesTree.ADD[d1][d2]) < epsilon)
		return 0;
	return 1;
}



//=================================================================
//==
//=================================================================
void copyInputTree(struct InputTree *tmpTree, struct InputTree aTree,int all,int allW){
	
	int n= aTree.size;
	int i,j;

	if(tmpTree->ADD == NULL){
		allocMemmory(tmpTree,n);
	}

	for(i=1;i<=2*n-1;i++)
		tmpTree->degre[i] = aTree.degre[i];

	for(i=1;i<=2*n-2;i++){
		for(j=1;j<=2*n-2;j++){
			tmpTree->Adjacence[i][j] = aTree.Adjacence[i][j];
			tmpTree->ADD[i][j] = aTree.ADD[i][j];
			if(all==1)
				tmpTree->Input[i][j] = aTree.Input[i][j];
		}
	}
  if(all==1){
  	 for(i=1;i<=2*n-2;i++){
		  for(j=1;j<=2*n-2;j++){
				tmpTree->Input[i][j] = aTree.Input[i][j];
		  }
	   }
	   
    for(i=1;i<=n;i++){	
			strcpy(tmpTree->SpeciesName[i],(const char*)aTree.SpeciesName[i]);
		}
	}
	
	if(allW==1){
    for(i=1;i<=n;i++){
		  for(j=1;j<=n;j++){
			 tmpTree->W[i][j] = aTree.W[i][j];
		  }
	 }
  }
	for(i=1;i<=2*n-3;i++){
		tmpTree->ARETE[2*i-1]  = aTree.ARETE[2*i-1];
		tmpTree->ARETE[2*i-2]  = aTree.ARETE[2*i-2];
		tmpTree->LONGUEUR[i-1] = aTree.LONGUEUR[i-1];
	}
	
	tmpTree->Root = aTree.Root;
	tmpTree->size = aTree.size;
	tmpTree->kt = aTree.kt;
}

void printLeaves(char *fichier,int nbHGT,struct HGT *outHGT,int noTree,char **NomsSpecies){

	FILE *out;
	int i,j;

	out = fopen(fichier,"a+");
	for(i=1;i<=nbHGT;i++){
		
		if((outHGT[i].valide == 1) || (outHGT[i].valide == 6)) {
			fprintf(out,"\ndetection %d : ",noTree);
			for(j=1;j<=outHGT[i].listSource[0];j++){
				fprintf(out,"%s",NomsSpecies[outHGT[i].listSource[j]]);
				if(j<outHGT[i].listSource[0])
					fprintf(out,";");
				else
					fprintf(out," ");
			}
			for(j=1;j<=outHGT[i].listDestination[0];j++){
				fprintf(out,"%s",NomsSpecies[outHGT[i].listDestination[j]]);
				if(j<outHGT[i].listDestination[0])
					fprintf(out,";");
			}
		}
	}
	fprintf(out,"\n");
	fclose(out);
}
int sameHGT(struct HGT HGT1, struct HGT HGT2){
	int i;

	if(HGT1.listSource[0] == HGT2.listSource[0] && HGT1.listDestination[0] == HGT2.listDestination[0]){
		
		for(i=1;i<=HGT1.listSource[0];i++){
			if(HGT1.listSource[i] != HGT2.listSource[i])
				return 0;
		}
		for(i=1;i<=HGT1.listDestination[0];i++){
			if(HGT1.listDestination[i] != HGT2.listDestination[i])
				return 0;
		}
		return 1;
	}

	return 0;
}

int sameHGT2(struct HGT HGT1, struct HGT HGT2){
	int i,h1=0,h2=0;

	if(HGT1.listSource[0] == HGT2.listSource[0]){ 
		h1 = 1;
		for(i=1;i<=HGT1.listSource[0];i++){
			if(HGT1.listSource[i] != HGT2.listSource[i])
				h1=0;
		}
	}

	if(HGT1.listDestination[0] == HGT2.listDestination[0]){
		h2=1;
		for(i=1;i<=HGT1.listDestination[0];i++){
			if(HGT1.listDestination[i] != HGT2.listDestination[i])
				h2 = 0;
		}
	}

	if(h2 == 1 && h1 == 1) return 1;


	if(HGT1.listSource[0] == HGT2.listDestination[0]){ 
		h1=1;
		for(i=1;i<=HGT1.listSource[0];i++){
			if(HGT1.listSource[i] != HGT2.listDestination[i])
				h1=0;
		}
	}

	if(HGT1.listDestination[0] == HGT2.listSource[0]){
		h2=1;
		for(i=1;i<=HGT1.listDestination[0];i++){
			if(HGT1.listDestination[i] != HGT2.listSource[i])
				h2 = 0;
		}
	}

	if(h2 == 1 && h1 == 1) return 1;
	
	return 0;
}

void copyHGT(struct HGT source,struct HGT *dest){
	int j;
	dest->trivial = source.trivial;
	dest->valide = source.valide;
	dest->source_A = source.source_A;
	dest->source_B = source.source_B;
	dest->dest_A = source.dest_A;
	dest->dest_B = source.dest_B;
	dest->crit.LS = source.crit.LS;
	dest->crit.RF = source.crit.RF;
	dest->crit.BD = source.crit.BD;
	dest->crit.rLS = source.crit.rLS;
	dest->crit.rRF = source.crit.rRF;
	dest->crit.rBD = source.crit.rBD;
	dest->crit.QD = source.crit.QD;
	dest->listSource = (int*) malloc((source.listSource[0]+1)*sizeof(int));
	dest->listDestination = (int*) malloc((source.listDestination[0]+1)*sizeof(int));
	dest->listSource[0] = source.listSource[0];
	for(j=1;j<=dest->listSource[0];j++){
		dest->listSource[j] = source.listSource[j];
	}
	dest->listDestination[0] = source.listDestination[0];
	for(j=1;j<=dest->listDestination[0];j++){
		dest->listDestination[j] = source.listDestination[j];
	}

}
//=============================================================================================================
//
//=============================================================================================================
void updateBootHGT(int first,struct HGT *bestHGT,int cpt_hgt,struct HGT *bootHGT, int *nbHGT_boot, double *tabBoot){
	
	int i,cpt=0,j;

	if(first==1){
		for(i=1;i<=cpt_hgt;i++){
			if((bestHGT[i].valide == 1) || (bestHGT[i].valide == 6)){
				cpt++;
				tabBoot[cpt] = 1;
				copyHGT(bestHGT[i],&bootHGT[cpt]);
			}
		}
		(*nbHGT_boot) = cpt;
	}
	else{
		for(i=1;i<=(*nbHGT_boot);i++){
			for(j=1;j<=cpt_hgt;j++){
				if((bestHGT[j].valide == 1) || (bestHGT[i].valide == 6))
					if (sameHGT(bootHGT[i],bestHGT[j])==1)
						tabBoot[i]+=1;
			}
		}
	}
	


}
//=====================================================
//
//=====================================================
void sortHGT(struct HGT *tabHGT,int nbHGT,struct Parameters param){
	
	int i,j;
	int criterion;
	int a,b,c,d,rf,qd;
	double ls,bd;

	if(strcmp(param.criterion,"rf")==0)
		criterion = 1;
	else if(strcmp(param.criterion,"ls")==0)
		criterion = 2;
	else if(strcmp(param.criterion,"bd")==0)
		criterion = 3;
	else
		criterion = 4;

	for(i=1;i<=nbHGT;i++){
		for(j=1;j<=nbHGT-i;j++){
			if( ((criterion == 1) && (tabHGT[j].crit.RF > tabHGT[j+1].crit.RF)) ||
				((criterion == 2) && (tabHGT[j].crit.LS > tabHGT[j+1].crit.LS)) ||
				((criterion == 4) && (tabHGT[j].crit.QD > tabHGT[j+1].crit.QD)) ||
				((criterion == 3) && (tabHGT[j].crit.BD > tabHGT[j+1].crit.BD)) ){
				rf = tabHGT[j].crit.RF;
				ls = tabHGT[j].crit.LS;
				bd = tabHGT[j].crit.BD;
				qd = tabHGT[j].crit.QD;
				a = tabHGT[j].source_A;
				b = tabHGT[j].source_B;
				c = tabHGT[j].dest_A;
				d = tabHGT[j].dest_B;
				tabHGT[j].crit.RF = tabHGT[j+1].crit.RF;
				tabHGT[j].crit.LS = tabHGT[j+1].crit.LS;
				tabHGT[j].crit.BD = tabHGT[j+1].crit.BD;
				tabHGT[j].crit.QD = tabHGT[j+1].crit.QD;
				tabHGT[j].source_A = tabHGT[j+1].source_A;
				tabHGT[j].source_B = tabHGT[j+1].source_B;
				tabHGT[j].dest_A = tabHGT[j+1].dest_A;
				tabHGT[j].dest_B = tabHGT[j+1].dest_B;
				tabHGT[j+1].crit.RF = rf;
				tabHGT[j+1].crit.LS = ls;
				tabHGT[j+1].crit.BD = bd;
				tabHGT[j+1].crit.QD = qd;
				tabHGT[j+1].source_A = a;
				tabHGT[j+1].source_B = b;
				tabHGT[j+1].dest_A = c;
				tabHGT[j+1].dest_B = d;
			}
		}
	}
}
//=========================================================
//
//=========================================================
void FreeMemory_InputTreeReduced(struct InputTree *aTree,int size){
	
	int i;
	int n=size;
	
	if(aTree->ADD != NULL){
	
		free(aTree->LONGUEUR);
		free(aTree->ARETE);

		for(i=0;i<2*aTree->size-1;i++){
			free(aTree->ADD[i]);
			free(aTree->W[i]);
		}
		for(i=0;i<2*n;i++)
			free(aTree->Adjacence[i]);

		if(aTree->SpeciesName != NULL){
			for(i=0;i<=n;i++){
				free(aTree->SpeciesName[i]);
			}
			free(aTree->SpeciesName);
		}
		
		free(aTree->ADD);
		free(aTree->Adjacence);
		free(aTree->W);
		free(aTree->degre);
	}
}

//=========================================================
//
//=========================================================
void FreeMemory_InputTree(struct InputTree *aTree,int size){
	
	int i;
	int n=size;
	
	if(aTree->ADD != NULL){
	
		free(aTree->LONGUEUR);
		free(aTree->ARETE);

		for(i=0;i<2*n;i++){
			free(aTree->ADD[i]);
			free(aTree->Input[i]);
			if(i<=n)
				free(aTree->W[i]);
		}
		for(i=0;i<2*n;i++){		
			free(aTree->Adjacence[i]);
		}

		if(aTree->SpeciesName != NULL){
			for(i=0;i<=n;i++){
				free(aTree->SpeciesName[i]);
			}
			free(aTree->SpeciesName);
		}
		
		free(aTree->ADD);
		free(aTree->Adjacence);
		free(aTree->Input);
		free(aTree->W);
		free(aTree->degre);
	}
}

//==============================
//
//==============================
int nextTreeIsNewick(FILE *in){
	
	char c;

	while(fscanf(in,"%c",&c)!= EOF){
		if(c != ' ' && c !='\t' && c != '\n' && c !='\r'){
			fseek(in,-1,SEEK_CUR);
			if(c=='(') return 1;
			return 0;
		}
	}
	
	return -1;
}

//==============================
//
//==============================
char * readNewick(FILE *in){
	
	int cpt=0;
	char c;
	char * newick = (char*)malloc(100000);

	do{
		c=(char)fgetc(in);
		newick[cpt] = c;
		cpt++;
	}while(c!=';');
	
	newick[cpt] = '\0';

	return newick;
}

//=======================================
//
//=======================================
int nbSpeciesNewick(string newick){
	int i=0;
	int n = 0;
	char symbol,parent=' ';
	char symbolOld =' ';  	
	int temoin =0;

	do{
		symbol = newick.at(i);
		i++;

		if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
}


//==================================================================
//
//==================================================================
void newickToMatrix(string newick,struct InputTree *aTree){
	
	int i,j;
	int pasBinaire=0;
	int pos_racine=-1;


	aTree->size = nbSpeciesNewick(newick);
	allocMemmory(aTree,aTree->size+1);
	pos_racine = lectureNewick(newick,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,&aTree->kt);

  {
    loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size,aTree->kt);
    Floyd(aTree->Adjacence,aTree->ADD,aTree->Input,aTree->size,aTree->kt); //1ere fois
  }
  
}

//===============================================
//
//===============================================
void readMatrix(FILE *in,struct InputTree *aTree){
	
	int n,i,j;
	char name[50];
	double value;

	fscanf(in,"%d",&n);
	aTree->size = n;
	allocMemmory(aTree,aTree->size+1);

	for(i=1;i<=n;i++){
		fscanf(in,"%s",name); strcpy(aTree->SpeciesName[i],name);
		for(j=1;j<=n;j++){
			fscanf(in,"%lf",&value);
			aTree->ADD[i][j] = aTree->Input[i][j] = value;
		}
	}
}

//===================================================================================================================================
//
//===================================================================================================================================
int readInputFile(string tree1, string tree2, const char *tmpFile, struct InputTree *speciesTree_t, struct InputTree *geneTree_t,char *fichier_erreur){
	string newick;
/* 	struct InputTree speciesTree_t;
	struct InputTree geneTree_t; */
	int ret;
	int finalTaille=0;
	
/* 	initInputTree(&speciesTree_t);
	initInputTree(&geneTree_t); */
	
	newick = tree1;
	/* cout<<"tree1 "<<newick<<endl; */
	newickToMatrix(newick,speciesTree_t);

	newick = tree2;
	/* cout<<"tree2 "<<newick<<endl; */
	newickToMatrix(newick,geneTree_t);
	
 	filtrerMatrice(speciesTree_t->Input,geneTree_t->Input,speciesTree_t->SpeciesName,geneTree_t->SpeciesName,speciesTree_t->size,geneTree_t->size,fichier_erreur);
	
	if((finalTaille=ecrireMatrice(speciesTree_t->Input,tmpFile,speciesTree_t->size,speciesTree_t->SpeciesName)) == -1)
		return -2;
	ajouterMatriceGene(geneTree_t->Input,tmpFile,geneTree_t->size,geneTree_t->SpeciesName);
	
	if(finalTaille<0){
		finalTaille = 0;
	}
	
	return finalTaille;
}




//=================================================================================
//
//=================================================================================
int midPoint(long int *ARETE,double **DIST,int n,int kt){

	double max;
	int i,j,i1,j1,P;


	/*/== recherche la plus grande distance*/
	max = 0;
	for(i=1;i<=n;i++){
		for(j=1;j<=n;j++){
			if(DIST[i][j] > max){
				max = DIST[i][j];
				i1=i; j1=j;
			}
		}
	}
	
	P = -1;
	for(i=1;i<=2*n-3-kt;i++){
		if(ARETE[2*i-1] == i1 || ARETE[2*i-1] == j1 || ARETE[2*i-2] == i1 || ARETE[2*i-2] == j1 ){
			if(DIST[ARETE[2*i-1]][ARETE[2*i-2]] >= max/2.0){
				P = i;
				break;
			}
		}
		else{
			if(ARETE[2*i-1] > n && ARETE[2*i-2]> n)
			if( ((DIST[i1][ARETE[2*i-1]] >= (max/2.0)) && (DIST[j1][ARETE[2*i-2]] >= (max/2.0)) ) ||
				((DIST[i1][ARETE[2*i-2]] >= (max/2.0)) && (DIST[j1][ARETE[2*i-1]] >= (max/2.0 ))) ){
					P = i;

					if( fabs(DIST[i1][j1] - DIST[i1][ARETE[2*i-1]] - DIST[ARETE[2*i-1]][j1]) < 0.00001 )
						break;
			}
		}
	}

	return P;
}

//=================================================================
//== 
//=================================================================
void AdjustBranchLength(struct InputTree *aTree1, struct InputTree aTree2,int binaire,int useFloyd){

	int i,j,iternumber=100;
	int kt = aTree1->kt;
	
	approx_arb(aTree2.ADD,aTree1->ADD,aTree1->ADD,aTree1->W,&iternumber,aTree1->ARETE,aTree1->LONGUEUR,1,&kt,binaire,aTree1->size);
  loadAdjacenceMatrix(aTree1->Adjacence,aTree1->ARETE, aTree1->LONGUEUR,aTree1->size,aTree1->kt);

	if(useFloyd == 1){
		Floyd(aTree1->Adjacence,aTree1->ADD,aTree1->size,aTree1->kt); // 3eme fois
  }
   
}

//=================================================================
//==
//=================================================================
void InitCriteria(struct CRITERIA * oldCrit, int size){

	int i;

	oldCrit->PLACE=(int *) malloc((2*size-3+1)*sizeof(int));
	oldCrit->PLACEI=(int *) malloc((2*size-3+1)*sizeof(int));
	oldCrit->B=(int **) malloc((2*size-3+1)*sizeof(int*));
	oldCrit->BI=(int **) malloc((2*size-3+1)*sizeof(int*));

	for (i=0;i<=2*size-3;i++)
	{
		oldCrit->B[i]=(int *) malloc((size+1)*sizeof(int));
		oldCrit->BI[i]=(int *) malloc((size+1)*sizeof(int));
	}

	oldCrit->BD=0;
	oldCrit->LS=0.0;
	oldCrit->RF=0;
}

void FreeCriteria(struct CRITERIA * Crit,int size){

	int i;

	free(Crit->PLACE);
	free(Crit->PLACEI);
	
	for(i=0;i<=2*size-3;i++){
		free(Crit->B[i]);
		free(Crit->BI[i]);
	}
	free(Crit->B);
	free(Crit->BI);
}

//====================================
//
//====================================

//================================================================================================================
//=
//================================================================================================================
int findAllMinimalScenario(struct InputTree SpeciesTree , struct InputTree GeneTree,int binaireSpecies,int binaireGene){

	return 0;
}

//================================================================================================================
//=
//================================================================================================================
int findAllHGT_no_criterion(struct InputTree SpeciesTree, struct InputTree GeneTree,struct Parameters param,struct HGT *tabHGT){
	
	int nbHGT=0;	
	struct InputTree tmpTree;
	int i,j,ktSpecies;
	int size = SpeciesTree.size;
	struct CRITERIA aCrit;
	struct DescTree *DTSpecies,*DTGene;

	initInputTree(&tmpTree);

	//== compute criteria
	InitCriteria(&aCrit,size);

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);

	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		for(j=i+1;j<2*size-3-SpeciesTree.kt;j++){
		
			//== is it a valid hgt ?
			if(isAValidHGT(SpeciesTree,i,j)==1){
				
				nbHGT++;
				tabHGT[nbHGT].source = i;
				tabHGT[nbHGT].destination = j;
				tabHGT[nbHGT].valide = 1;

			}
		}
	
	FreeCriteria(&aCrit,size);

	for(i=1;i<=nbHGT;i++)
		findListSpecies(&tabHGT[i],DTSpecies,SpeciesTree);


	return nbHGT;
}

//================================================================================================================
//=
//================================================================================================================
int findAllHGT(struct InputTree SpeciesTree, struct InputTree GeneTree,struct Parameters param,struct HGT *tabHGT){
	
	int nbHGT=0;	
	struct InputTree tmpTree;
	int i,j,ktSpecies;
	int size = SpeciesTree.size;
	struct CRITERIA aCrit;
	struct DescTree *DTSpecies,*DTGene;

	initInputTree(&tmpTree);

	//== compute criteria
	InitCriteria(&aCrit,size);

	DTGene = (struct DescTree*)malloc((2*GeneTree.size)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size)*sizeof(struct DescTree));

	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		for(j=1;j<2*size-3-SpeciesTree.kt;j++){
		
			//== is it a valid hgt ?
			if(isAValidHGT(SpeciesTree,i,j)==1){

				copyInputTree(&tmpTree,SpeciesTree,0,0);

				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,i,j,DTSpecies,DTGene) == 0) continue;
				
				applyHGT(GeneTree.ADD,&tmpTree,i,j);
				AdjustBranchLength(&tmpTree,GeneTree,0,1);
				
				computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				nbHGT++;
				loadCriteria(aCrit,&tabHGT[nbHGT]);
				tabHGT[nbHGT].source = i;
				tabHGT[nbHGT].destination = j;
				
			}
		}
	
	FreeCriteria(&aCrit,size);

	for(i=1;i<=nbHGT;i++)
		findListSpecies(&tabHGT[i],DTSpecies,SpeciesTree);

	return nbHGT;
}
//===============================================================================================================
//== 
//===============================================================================================================
int findBestHGT(int initial,struct InputTree SpeciesTree,struct InputTree GeneTree,struct Parameters param,struct HGT *aHGT){

	struct InputTree tmpTree;
	int i,j,k,l,first=1,ret=0;
	int size = SpeciesTree.size;
	int ktSpecies;
	struct CRITERIA aCrit,aCritRef;
	struct DescTree *DTSpecies,*DTGene;

	initInputTree(&tmpTree);

	//== compute criteria
	InitCriteria(&aCrit,size);
	InitCriteria(&aCritRef,size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	loadCriteria(aCrit,aHGT);

	DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
	
	//== Test all possible HGT
	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		for(j=1;j<2*size-3-SpeciesTree.kt;j++){
			
			//== is it a valid hgt ?
			if(isAValidHGT(SpeciesTree,i,j)==1 && i!=j ){

				copyInputTree(&tmpTree,SpeciesTree,0,0);
				
				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,i,j,DTSpecies,DTGene) == 0) {
						continue;
						tmpTree.ADD = NULL;
						tmpTree.ARETE = NULL;
						
					}
				
					applyHGT(GeneTree.ADD,&tmpTree,i,j);
		
					AdjustBranchLength(&tmpTree,GeneTree,0,1);
					if(initial == 1){
						computeCriteria(tmpTree.ADD,SpeciesTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,SpeciesTree.LONGUEUR,SpeciesTree.ARETE);
						if(aCrit.RF == 1){
							ret += TestCriterionAndUpdate(&first,param.criterion,aCrit,aHGT,i,j,0,0);
							i=j=(int)INFINI;
							break;
						}
					}else{
						computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
						ret += TestCriterionAndUpdate(&first,param.criterion,aCrit,aHGT,i,j,0,0);
						if(aHGT->crit.RF == 0 /*|| aHGT->crit.LS < epsilon */|| aHGT->crit.BD == 0 || aHGT->crit.QD == 0) { i=j=(int)INFINI;}
					}

			}
		}
	if(ret > 0 )
		findListSpecies(aHGT,DTSpecies,SpeciesTree);
	
	deleteBipartition(DTSpecies,SpeciesTree);
	deleteBipartition(DTGene,GeneTree);
	FreeCriteria(&aCrit,size);
	FreeCriteria(&aCritRef,size);
	FreeMemory_InputTree(&tmpTree,tmpTree.size);

	return ret;
}

//===========================================================================================================================
//== 
//===========================================================================================================================
int findBestHGT_nombreLimite(struct DescTree *DTSpecies,struct DescTree *DTGene,int * tab_branches,int nb_branches,struct InputTree SpeciesTree,struct InputTree GeneTree,struct Parameters param,struct HGT *aHGT){

	struct InputTree tmpTree;
	int i,j,k,l,first=1,ret=0;
	int size = SpeciesTree.size;
	int ktSpecies;
	struct CRITERIA aCrit,aCritRef;
	
	initInputTree(&tmpTree);

	//== compute criteria
	InitCriteria(&aCrit,size);
	InitCriteria(&aCritRef,size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	loadCriteria(aCrit,aHGT);

	//== Test all possible HGT
	for(i=1;i<=nb_branches;i++){
		for(j=1;j<=nb_branches;j++){
			if(isAValidHGT(SpeciesTree,tab_branches[i],tab_branches[j])==1 && i!=j ){
			
				copyInputTree(&tmpTree,SpeciesTree,0,0);

				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,tab_branches[i],tab_branches[j],DTSpecies,DTGene) == 0) {
						continue;
						
						tmpTree.ADD = NULL;
						tmpTree.ARETE = NULL;
						
					}

				applyHGT(GeneTree.ADD,&tmpTree,tab_branches[i],tab_branches[j]);
				AdjustBranchLength(&tmpTree,GeneTree,0,1);
	
				computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);

				ret += TestCriterionAndUpdate(&first,param.criterion,aCrit,aHGT,tab_branches[i],tab_branches[j],0,0);
				if(aHGT->crit.RF == 0 /*|| aHGT->crit.LS < epsilon */|| aHGT->crit.BD == 0 || aHGT->crit.QD == 0) { i=j=(int)INFINI;}

			}
		}
		
	}
	
	FreeCriteria(&aCrit,size);
	FreeCriteria(&aCritRef,size);
	FreeMemory_InputTree(&tmpTree,tmpTree.size);

	return ret;
}


//===============================================================================================================
//== 
//===============================================================================================================
int findBestHGTtab(struct InputTree SpeciesTree,struct InputTree GeneTree,struct Parameters param,struct HGT *aHGT, int *nbHgtFound,int *initial,int bootstrap){

	struct InputTree tmpTree;
	struct InputTree tmpTree2;
	
	int i,j,k,l,first=1,ret=0,trouve;
	int size = SpeciesTree.size;
	int ktSpecies;
	struct CRITERIA aCrit,aCritRef,aCritRef2,aCrit_tmpTree2;
	struct DescTree *DTSpecies,*DTGene;
	int encore = 0;
  
	initInputTree(&tmpTree);
	initInputTree(&tmpTree2);
	
	(*nbHgtFound) = 0;

	printf("\n\n== NEW STEP OF DETECTION == [size of the trees (after reduction)=%d]",SpeciesTree.size);
  
	//== compute criteria
	InitCriteria(&aCrit,size);
	InitCriteria(&aCritRef,size);
	InitCriteria(&aCritRef2,size);
	InitCriteria(&aCrit_tmpTree2,size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef2,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	loadCriteria(aCrit,&aHGT[0]);
  
	DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);

	int flag2 = 1;
	int flag3;   		//== pour tester les especes avec sur la branche au dessus
	int initial3=1;
	
	
	do{
		for(i=1;i<2*size-3-SpeciesTree.kt;i++)
			for(j=i+1;j<2*size-3-SpeciesTree.kt;j++){
				//== is it a valid hgt ?
				trouve=0;
				if(isAValidHGT(SpeciesTree,i,j)==1 && i!=j ){
					copyInputTree(&tmpTree,SpeciesTree,0,0);
				
					if(strcmp(param.subtree,"yes") == 0){
						if(TestSubTreeConstraint(SpeciesTree,i,j,DTSpecies,DTGene) == 0)  {
							continue;
							tmpTree.ADD = NULL;
							tmpTree.ARETE = NULL;
						}
						flag3=0;
						if(TestSubTreeLeafs(SpeciesTree,i,j,DTSpecies,DTGene) == 1){
							flag3=1;
						}
					}
					
					applyHGT(GeneTree.ADD,&tmpTree,i,j);
					
					if(strcmp(param.criterion,"ls") == 0)
						AdjustBranchLength(&tmpTree,GeneTree,0,1);
					computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				
					first=1;
					loadCriteria(aCritRef,&aHGT[(*nbHgtFound)]);
					if(((aCritRef.RF-aCrit.RF) == 1) && (*initial==1) && (
					  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-1]) || 
					  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-2]) || 
					  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2]) || 
					  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2])
					  )) {
						UpdateCriterion(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],i,j,0);
						
						aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*i-1];
						aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*i-2];
						aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*j-1];
						aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*j-2];
						
						findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
						trouve=1;
					}
					else  if((*initial==0) && ((flag3 == 1) || (initial3 == 0) )){
						if(TestCriterionAndUpdate(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],i,j,0,0) == 1){		
							aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*i-1];
							aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*i-2];
							aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*j-1];
							aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*j-2];
							findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
							trouve=1;
						}
					}
				}

				if(isAValidHGT(SpeciesTree,j,i)==1 && i!=j ){
					copyInputTree(&tmpTree2,tmpTree,0,0);
					copyInputTree(&tmpTree,SpeciesTree,0,0);
					
					if(strcmp(param.subtree,"yes") == 0){
						if(TestSubTreeConstraint(SpeciesTree,j,i,DTSpecies,DTGene) == 0) {
							continue;
							tmpTree.ADD = NULL;
							tmpTree.ARETE = NULL;
						}
					
						if(TestSubTreeLeafs(SpeciesTree,j,i,DTSpecies,DTGene)==1){
							flag3=2;
						 }
					}
					
					applyHGT(GeneTree.ADD,&tmpTree,j,i);
					
					if(strcmp(param.criterion,"ls") == 0)
					  AdjustBranchLength(&tmpTree,GeneTree,0,1);
					computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
					if(tmpTree2.ADD == NULL){
						printf("HOUSTON, we have a problem !");
						aCrit_tmpTree2.RF = (int) INFINI;
					}
					else{
						computeCriteria(tmpTree.ADD,tmpTree2.ADD,size,&aCrit_tmpTree2,tmpTree.LONGUEUR,tmpTree.ARETE,tmpTree2.LONGUEUR,tmpTree2.ARETE);
					}
					int flag=0;
					first=1;
					int flag_bootstrap = 0;
					if(trouve==0)
						loadCriteria(aCritRef,&aHGT[(*nbHgtFound)]);
					else{
						//printf("\n===========> rand_bootstrap : RF = %d",aCrit_tmpTree2.RF);
						if((bootstrap == 1) && (aCrit_tmpTree2.RF == 0) && /*(aCrit.BD == aHGT[(*nbHgtFound)].crit.BD) && */(aCrit.RF == aHGT[(*nbHgtFound)].crit.RF) ){
							flag_bootstrap = 1;
							//printf("\n==> rand_bootstrap1 : [%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)].source_A,aHGT[(*nbHgtFound)].source_B,aHGT[(*nbHgtFound)].dest_A,aHGT[(*nbHgtFound)].dest_B);
						}
						else{
							aHGT[(*nbHgtFound)].crit.diff_bd = fabs(aCrit.BD - aHGT[(*nbHgtFound)].crit.BD);
							if((flag2 == 1) && (fabs(aCrit.BD - aHGT[(*nbHgtFound)].crit.BD) <= 1)){
								if( (aHGT[(*nbHgtFound)].crit.RF - aCrit.RF) >= 1 ){
									trouve = 0;
									flag = 2;
								}      
								else{
									flag = 1;
									//printf("\n--> TEST FLAG : RF1=%d | RF2=%d | BD1=%2.1lf | BD2=%2.1lf",aHGT[(*nbHgtFound)].crit.RF,aCrit.RF,aHGT[(*nbHgtFound)].crit.BD,aCrit.BD);
								}
							}
							else{
								//printf("\n-->[%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)].source_A,aHGT[(*nbHgtFound)].source_B,aHGT[(*nbHgtFound)].dest_A,aHGT[(*nbHgtFound)].dest_B);
								//printf("\tRF1=%d | RF2=%d | BD1=%2.1lf | BD2=%2.1lf",aHGT[(*nbHgtFound)].crit.RF,aCrit.RF,aHGT[(*nbHgtFound)].crit.BD,aCrit.BD);
							}
						}
					}
					
					//flag_bootstrap = 0;
					if(flag != 1){
						if(((aCritRef.RF-aCrit.RF) == 1) && (*initial==1) && (
						  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-1]) || 
						  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-2]) || 
						  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2]) || 
						  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2])
						)) {
						
							if ((flag_bootstrap == 1) && (rand_bootstrap == 1)){
								
								//printf("\n===========>REVERSE : initial = %d",*initial);
								//printf("\tRF1=%d | RF2=%d | BD1=%2.1lf | BD2=%2.1lf",aHGT[(*nbHgtFound)].crit.RF,aCrit.RF,aHGT[(*nbHgtFound)].crit.BD,aCrit.BD);
								UpdateCriterion(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],j,i,flag_bootstrap);
								aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*j-1];
								aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*j-2];
								aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*i-1];
								aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*i-2];
								findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
								trouve=1;
								//printf("\n==> rand_bootstrap1.5 : [%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)].source_A,aHGT[(*nbHgtFound)].source_B,aHGT[(*nbHgtFound)].dest_A,aHGT[(*nbHgtFound)].dest_B);
							}
						}
						else if(((*initial==0) && ((flag3 == 2) || (initial3 == 0))) || (flag_bootstrap == 1)){
							if(TestCriterionAndUpdate(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],j,i,flag,flag_bootstrap) == 1){	
								//printf("\n===========>REVERSE : flag3 = %d, initial3=%d",flag3,initial3);
								//printf("\tRF1=%d | RF2=%d | BD1=%2.1lf | BD2=%2.1lf",aHGT[(*nbHgtFound)].crit.RF,aCrit.RF,aHGT[(*nbHgtFound)].crit.BD,aCrit.BD);
								aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*j-1];
								aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*j-2];
								aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*i-1];
								aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*i-2];
								findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
								trouve=1;
								//printf("\n==> rand_bootstrap2 : [%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)].source_A,aHGT[(*nbHgtFound)].source_B,aHGT[(*nbHgtFound)].dest_A,aHGT[(*nbHgtFound)].dest_B);
							}
						}
					}				
				}
				
				if (trouve == 1){
					(*nbHgtFound)++;
					//printf("\n===================> rand_bootstrap final : [%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)-1].source_A,aHGT[(*nbHgtFound)-1].source_B,aHGT[(*nbHgtFound)-1].dest_A,aHGT[(*nbHgtFound)-1].dest_B);			
				}

			}
			encore = 0;
		
		if((*nbHgtFound == 0) && (flag2==1)){
		  flag2 = 0;
		  encore = 1;
		}
		
		if((*nbHgtFound == 0) && (initial3==1) && (*initial==0) ){
		  initial3 = 0;
		  flag2=1;
		  encore = 1;
		}
		
		if((*nbHgtFound == 0) && (*initial==1)){
		  (*initial) = 0;
		  encore = 1;
		  flag2=1;
		}
		
		//initial3 = 0;

		
	 //   printf("On est rendu icit");
	}while(encore == 1);
	
	
	
	deleteBipartition(DTSpecies,SpeciesTree);
	deleteBipartition(DTGene,GeneTree);
	FreeCriteria(&aCrit,size);
	FreeCriteria(&aCritRef,size);
	FreeMemory_InputTree(&tmpTree,tmpTree.size);

	//printf("\n3)RF=%d",aHGT[0].crit.RF);	
	return (*nbHgtFound);
}



//=================================================================
//==
//=================================================================
void findBranch(struct InputTree aTree,int *branch, int * elt){
	
	double max = INFINI; //= disctance between root and intersection
	int j;

	(*branch) = 0;
	
	if(elt[0]==1){
		for(j=1;j<=2*aTree.size-3-aTree.kt;j++)
			if(aTree.ARETE[2*j-1]==elt[1] || aTree.ARETE[2*j-2]==elt[1])
				(*branch) = j;
	}
	else{
		for(j=2;j<=elt[0];j++){		
			if(max > ( (aTree.ADD[aTree.Root][elt[1]] + aTree.ADD[aTree.Root][elt[j]] - aTree.ADD[elt[1]][elt[j]]) / 2.0 ) )
				max = (aTree.ADD[aTree.Root][elt[1]] + aTree.ADD[aTree.Root][elt[j]] - aTree.ADD[elt[1]][elt[j]]) / 2.0 ;  
		}
		//printf("\nmax=%lf",max);	

		for(j=1;j<=2*aTree.size-3-aTree.kt;j++){
			if( (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]]- max) < 2*epsilon) && (aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]] < aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]]) && (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]]+aTree.ADD[aTree.ARETE[2*j-1]][elt[1]] - aTree.ADD[aTree.Root][elt[1]]) < 2*epsilon)  )
				(*branch) = j;
			if( (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]]- max) < 2*epsilon) && (aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]] < aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]]) && (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]]+aTree.ADD[aTree.ARETE[2*j-2]][elt[1]] - aTree.ADD[aTree.Root][elt[1]]) < 2*epsilon))
				(*branch) = j;
		}
	}
}

//======================================================
//==
//======================================================
void expandBestHGT(struct HGT bestHGTRed,struct HGT *bestHGT,struct ReduceTrace aMap,struct DescTree * DTSpecies,struct InputTree SpeciesTree){
	
	//bestHGT->source = (int*)malloc(

	int i,j,nbSource=0,nbDest=0,tmp;

	for(i=1;i<=bestHGTRed.listSource[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listSource[i]]].nbSommet;j++){
			nbSource++;
		}
	}

	bestHGT->listSource = (int*)malloc((nbSource+1)*sizeof(int));
	bestHGT->listSource[0] = nbSource;
	tmp=1;
	for(i=1;i<=bestHGTRed.listSource[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listSource[i]]].nbSommet;j++){
			bestHGT->listSource[tmp] = DTSpecies[aMap.map[bestHGTRed.listSource[i]]].Tableau[j];
			tmp++;
		}
	}
	TrierTableau(bestHGT->listSource,tmp-1);
	
/*	printf("\nSource :");
	for(i=1;i<=bestHGT->listSource[0];i++)
		printf(" %d",bestHGT->listSource[i]);
    */

	for(i=1;i<=bestHGTRed.listDestination[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listDestination[i]]].nbSommet;j++){
			nbDest++;
		}
	}

	bestHGT->listDestination = (int*)malloc((nbDest+1)*sizeof(int));
	bestHGT->listDestination[0] = nbDest;
	tmp=1;
	for(i=1;i<=bestHGTRed.listDestination[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listDestination[i]]].nbSommet;j++){
			bestHGT->listDestination[tmp] = DTSpecies[aMap.map[bestHGTRed.listDestination[i]]].Tableau[j];
			tmp++;
		}
	}
	TrierTableau(bestHGT->listDestination,tmp-1);
	
/*	printf("\nDest   :");
	for(i=1;i<=bestHGT->listDestination[0];i++)
		printf(" %d",bestHGT->listDestination[i]);
	*/

	//== find branches
	findBranch(SpeciesTree,&i,bestHGT->listSource);
	findBranch(SpeciesTree,&j,bestHGT->listDestination);
	bestHGT->source = i;
	bestHGT->destination = j;

	bestHGT->source_A = SpeciesTree.ARETE[2*i-1];
	bestHGT->source_B = SpeciesTree.ARETE[2*i-2];
	bestHGT->dest_A = SpeciesTree.ARETE[2*j-1];
	bestHGT->dest_B = SpeciesTree.ARETE[2*j-2];

	bestHGT->valide = bestHGTRed.valide;
	bestHGT->crit = bestHGTRed.crit;

	free(bestHGTRed.listSource);
	free(bestHGTRed.listDestination);
	bestHGTRed.listSource = NULL;
	bestHGTRed.listDestination = NULL;
}

//======================================================
//==
//======================================================
int isInside(struct DescTree DT1,struct DescTree DT2){

	int i;

	for(i=1;i<=DT1.nbSommet;i++){
		if(DT2.Tableau[1] == DT1.Tableau[i]) return 1;
	}
	return 0;
}

//=================================================================
//== 
//=================================================================
void CreateSubStructures(struct InputTree * aTree,int inc,int binaire){
	
	int n = aTree->size;
	int i,j;
	int kt=0;
	inc = 10;
		
	//printf("n=%d",n);
	if(aTree->ARETE == NULL){
		aTree->ARETE    =(long int*)malloc(4*(2*(n+inc))*sizeof(long int));
		aTree->LONGUEUR	=(double*)malloc((4*(n+inc))*sizeof(double));
		aTree->Adjacence=(double**)malloc((2*(n+inc)+1)*sizeof(double*));
		for(i=0;i<2*(n+inc);i++)
			aTree->Adjacence[i]=(double*)malloc((2*(n+inc)+1)*sizeof(double));
	}
//	printf("\nbinaire = %d",binaire);
	kt = aTree->kt = Tree_edges (aTree->ADD,aTree->ARETE,aTree->LONGUEUR,n,binaire);
	
	//printf("\n\n");
	//printEdges(0,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,aTree->size);
	//for(i=1;i<=2*n-3-kt;i++)
	//	printf("\n%d-%d : %lf",aTree->ARETE[2*i-1],aTree->ARETE[2*i-2],aTree->LONGUEUR[i-1]);
	//scanf("%d",&i);
	loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,n,aTree->kt);
	Floyd(aTree->Adjacence,aTree->ADD,n,aTree->kt); // 4eme fois
//  global_cpt4++;
	
  //===creation de degre
	aTree->degre = (int*)malloc(2*(n+inc)*sizeof(int));
	for(i=1;i<=2*n-2-kt;i++){
		aTree->degre[i]=0;
		for(j=1;j<=2*n-2-kt;j++)
			if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
	}
	/*for(i=1;i<=2*n-2-kt;i++)
		printf("%d ",aTree->degre[i]);*/
	
}

//=================================================================
//==
//=================================================================
void ReduceTree(struct InputTree SpeciesTree,struct InputTree GeneTree,struct InputTree *SpeciesTreeRed,struct InputTree *GeneTreeRed,struct ReduceTrace *aTrace,struct DescTree *DTSpecies,struct DescTree *DTGene,int binaireSpecies,int binaireGene){


	int i,j,l=1,tmp,pos; //k;
	struct CRITERIA aCrit;

	InitCriteria(&aCrit,SpeciesTree.size);

	aTrace->species  = (int *)malloc((2*SpeciesTree.size)*sizeof(int));
	aTrace->gene     = (int *)malloc((2*SpeciesTree.size)*sizeof(int));
	aTrace->map		 = (int *)malloc((2*SpeciesTree.size)*sizeof(int));
	//aTrace->speciesToGene = (int *)malloc((2*SpeciesTree.size-1)*sizeof(int));
	aTrace->size = 0;

	for(i=0;i<2*SpeciesTree.size;i++){
		aTrace->species[i] = 0;
		aTrace->gene[i] = 0;
	}
	
	for(i=1;i<2*(SpeciesTree.size)-2-SpeciesTree.kt;i++)
		for(j=1;j<2*(GeneTree.size)-2;j++){
		
			if(DTSpecies[i].nbSommet == DTGene[j].nbSommet && DTSpecies[i].nbSommet > 1 && i!= SpeciesTree.size && j!= GeneTree.size){
				if(vecteursEgaux(DTSpecies[i],DTGene[j]) == 1){
				
						computeCriteria(DTSpecies[i].Matrice,DTGene[j].Matrice,DTSpecies[i].nbSommet+1,&aCrit,NULL,NULL,NULL,NULL);
						if(aCrit.RF == 0){
							aTrace->species[i] = 1;
							aTrace->gene[j] = 1;
						}
				}
			}
			if(DTSpecies[i].nbSommet == 1 ){
				aTrace->species[i] = 1;
				aTrace->gene[i] =1;
				//aTrace->speciesToGene[i] = i;
			}
		}

		
		//== ici on test les sous-ensembles 2 par 2, si un est le sous-ensemble de l'autre, il n'est pas necessaire de consid�rer
		//== ce dernier.
		for(i=1;i<=2*SpeciesTree.size-2-SpeciesTree.kt;i++){
			for(j=1;j<=2*SpeciesTree.size-2-SpeciesTree.kt;j++){
				if(aTrace->species[i] != 0 && aTrace->species[j] != 0 && i!=j && i!=SpeciesTree.Root && j!=SpeciesTree.Root)
					if(DTSpecies[i].nbSommet > DTSpecies[j].nbSommet && isInside(DTSpecies[i],DTSpecies[j])==1){
						aTrace->species[j] = 0;
					}
			}
		}



		//== creer les matrices de distances
		tmp=0;
		pos=0;
		for(i=1;i<=2*SpeciesTree.size-2-SpeciesTree.kt;i++){
			if(aTrace->species[i] != 0){
				pos++;
				aTrace->map[pos] = i; //DTSpecies[i].Tableau[1];
	//			printf("<br> %d-%d",pos,i);
			}
		}
	/*	printf("\n");
		for(i=1;i<=pos;i++){
			printf("\n%d-->%d (%d)",i,aTrace->map[i],DTSpecies[aTrace->map[i]].Tableau[1]);
		}*/

		pos++; //= racine
		if(SpeciesTreeRed->ADD == NULL){
			SpeciesTreeRed->ADD = (double**)malloc((2*pos-1)*sizeof(double*));
			GeneTreeRed->ADD = (double**)malloc((2*pos-1)*sizeof(double*));
			SpeciesTreeRed->W = (double**)malloc((2*pos-1)*sizeof(double*));
			GeneTreeRed->W = (double**)malloc((2*pos-1)*sizeof(double*));
			for(i=0;i<2*pos-1;i++){
				SpeciesTreeRed->W[i] = (double*)malloc((2*pos-1)*sizeof(double));
				GeneTreeRed->W[i] = (double*)malloc((2*pos-1)*sizeof(double));
			}
			for(i=0;i<2*pos-1;i++){
				SpeciesTreeRed->ADD[i] = (double*)malloc((2*pos-1)*sizeof(double));
				GeneTreeRed->ADD[i] = (double*)malloc((2*pos-1)*sizeof(double));
			}
		}
		for(i=1;i<pos;i++)
			if(aTrace->map[i] != SpeciesTree.Root){
				for(j=1;j<pos;j++){
					if(aTrace->map[i] != SpeciesTree.Root){
						SpeciesTreeRed->ADD[i][j] = SpeciesTree.ADD[DTSpecies[aTrace->map[i]].Tableau[1]][DTSpecies[aTrace->map[j]].Tableau[1]];
						GeneTreeRed->ADD[i][j] = GeneTree.ADD[DTSpecies[aTrace->map[i]].Tableau[1]][DTSpecies[aTrace->map[j]].Tableau[1]];
					}
				}
			}
		
		aTrace->species[SpeciesTree.Root] = 1;
		aTrace->map[pos] = SpeciesTree.Root;
		for(i=1;i<pos;i++){
			SpeciesTreeRed->ADD[i][pos] = SpeciesTreeRed->ADD[pos][i] = SpeciesTree.ADD[aTrace->map[pos]][DTSpecies[aTrace->map[i]].Tableau[1]];
			GeneTreeRed->ADD[i][pos] = GeneTreeRed->ADD[pos][i] = GeneTree.ADD[aTrace->map[pos]][DTSpecies[aTrace->map[i]].Tableau[1]];
		}
		
		SpeciesTreeRed->size = GeneTreeRed->size = pos;
		SpeciesTreeRed->Root = GeneTreeRed->Root = pos;
		
		CreateSubStructures(SpeciesTreeRed,0,binaireSpecies);
		CreateSubStructures(GeneTreeRed,0,binaireGene);

		FreeCriteria(&aCrit,SpeciesTree.size);
	//	printf("\n");
	//	printMatrix(SpeciesTree.SpeciesName,SpeciesTreeRed->ADD,SpeciesTreeRed->size);
}

//=================================================================
//==
//=================================================================
int DeleteUseLessHGT(int nbHGT,struct HGT *bestHGT, struct InputTree SpeciesTree, struct InputTree SaveTree){

	int i,j,source,destination,p,q;
	struct InputTree tmpTree;
	struct CRITERIA aCrit;
	int cpt=0,retour=0;

	InitCriteria(&aCrit,SpeciesTree.size);
	initInputTree(&tmpTree);
	printf("\nHGT_INT : Procedure de suppression des transferts inutiles, il y a %d transferts",nbHGT);
	for(i=1;i<=nbHGT;i++){
		printf("\nHGT_INT : bestHGT[%d].valide=%d",i,bestHGT[i].valide);
	}
	for(i=1;i<=nbHGT;i++){
		if((bestHGT[i].valide == 1) || (bestHGT[i].valide == 6) || (bestHGT[i].valide == 4))
		{
			copyInputTree(&tmpTree,SaveTree,0,0);
			for(j=1;j<=nbHGT;j++){
				if((bestHGT[j].valide == 1) || (bestHGT[j].valide == 6) || (bestHGT[j].valide == 4))
				{
					if (i!=j){
						findBranch(tmpTree,&source,bestHGT[j].listSource);
						findBranch(tmpTree,&destination,bestHGT[j].listDestination);
						if(source != 0 && destination !=0){
							if(isAValidHGT(tmpTree,source,destination)==1){
								applyHGT(SpeciesTree.ADD,&tmpTree,source,destination);
							}
						}
					}
				}
			}
			AdjustBranchLength(&tmpTree,SpeciesTree,0,1);
			//printf("\nHGT-DETECTION : on calcule les criteres");
			computeCriteria(SpeciesTree.ADD,tmpTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);
			printf("\nHGT-DETECTION : on a termine le calcul RF=%d",aCrit.RF);
			if(aCrit.RF == 0){
				printf("\nHGT_INT : Le transfert %d est inutile",i);
				bestHGT[i].valide = 0;
				retour++;
			}
		}
	}
	FreeMemory_InputTree(&tmpTree,tmpTree.size);
	FreeCriteria(&aCrit,SpeciesTree.size);
	return retour;
}




//=================================================================
//==
//=================================================================
void help(){
	printf("\nHGT-DETECTION version 3.0");
	printf("\nby Alix Boc and Vladimir Makarenkov");

	printf("\n\nUsage :\nhgt -inputfile=[inputfilename] -outputfile=[outputfilename] -criterion=[rf|ls|bd]");
	printf(" -version=[web|consol] -speciesroot=[midpoint|prompt|file] -generoot=[midpoint|prompt|file]");
	printf(" -load=[no|yes] -viewTree=[no|yes] -scenario=[unique|multiple] -nbhgt=[maxhgt] -path=[path]");

	printf("\n\ncriterion          [rf] = Robinson and Foulds distance (default)");
	printf("  \n                   [ls] = Least-Square optimization");
	printf("  \n                   [bd] = Bipartition distance");
	printf("\n\nversion            [consol] (default)");
	printf("  \n                   [web] = get result file and tree files in web ");
	printf("  \n                   format for the Trex-online web site");
	printf("\n\nspeciesroot        [midpoint] = the root is selected by the midpoint ");
	printf("  \n                   method (default)");
	printf("  \n                   [prompt] = the program ask for the root branch");
	printf("  \n                   [file] = the root is in a file called speciesroot.txt");
	printf("\n\ngeneroot           [midpoint] = the root is selected by the midpoint ");
	printf("  \n                   method (default)");
	printf("  \n                   [prompt] = the program ask for the root branch");
	printf("  \n                   [file] = the root is in a file called generoot.txt");
	printf("\n\nsubtree            [yes] = use the subtree constraint (default)");
	printf("  \n                   [no]");
	printf("\n\nscenario           [unique] (default)");
	printf("  \n                   [multiple]");
	printf("\n\nnbhgt              number max of hgts for unique and multiple scenario. ");
	printf("  \n                   Default value = 50 hgts");
	printf("\n\nload               [yes] = load the tree from files speciestree.txt ");
	printf("  \n                   genetree.txt speciesroot.txt generoot.txt");
	printf("  \n                   [no] (default)");
	printf("\n\npath               path without \"/\" at the end. Default path is \".\"");
	printf("  \n                   the path will be add to the file. ex : path/input.txt");

	printf("\n\nExample. \nCompute hgt-detection with default parameters : ./hgt -inputfile=input.txt\n\n");

}


void printTransfer(FILE *out,int sortie,char ** noms,int n,int source_A,int source_B,int dest_A,int dest_B){

	if(sortie == 0)
		fprintf(out,"\nFrom branch %d--%d to branch %d--%d",source_A,source_B,dest_A,dest_B);
	else{
		if(source_A	<= n)
			fprintf(out,"\nFrom branch %s--%d",noms[source_A],source_B);
		else if(source_B <= n)
			fprintf(out,"\nFrom branch %d--%s",source_A,noms[source_B]);
		else
			fprintf(out,"\nFrom branch %d--%d",source_A,source_B);

		if(dest_A <= n)
			fprintf(out," to branch %s--%d",noms[dest_A],dest_B);
		else if(dest_B <= n)
			fprintf(out," to branch %d--%s",dest_A,noms[dest_B]);
		else
			fprintf(out," to branch %d--%d",dest_A,dest_B);
	}
}

void printListSpecies(struct HGT * bestHGT, int nbHGT, double *tabBoot, int bootmin){

		FILE * out = fopen ("result.txt","w+");
		for(int i=1; i<= nbHGT && (bestHGT[i].valide == 1 || (bestHGT[i].valide == 6)); i++){
			if((tabBoot != NULL && tabBoot[i] >= bootmin) || (tabBoot == NULL)){
				fprintf(out,"\n%d",bestHGT[i].listSource[0]);
				for(int j=1;j<=bestHGT[i].listSource[0];j++)
					fprintf(out," %d",bestHGT[i].listSource[j]);
				fprintf(out,"\n%d",bestHGT[i].listDestination[0]);
				for(int j=1;j<=bestHGT[i].listDestination[0];j++)
					fprintf(out," %d",bestHGT[i].listDestination[j]);
			}
		}
		fclose(out);
}
//==============================================================================
//== AFFICHAGE DES RESULTATS DANS LE FICHIER DE SORTIE
//==============================================================================
void printHGT(FILE * res,FILE * res2,char *stepbystep,struct CRITERIA * multicheckTab,char *mode,int RFref,/*FILE *out,*/struct InputTree SpeciesTree,struct HGT *tabHGT, int nbHGT, double *tabBoot,char *subtree,int bootmin){
	int i,j,k;
	int diffHGT=0;
	FILE * res1 = res; //fopen ("result2.txt","w+");
	int nbTrivial=0;
	
	fprintf(res,"mode=%s\n",mode);
	
/*   if(strcmp(mode,"monocheck") == 0 || multicheckTab == NULL){
		for(i=1;i<=nbHGT;i++){
			if(i==1)
				diffHGT = RFref - tabHGT[i].crit.RF;
			else
				diffHGT = tabHGT[i-1].crit.RF - tabHGT[i].crit.RF;
		//	fprintf(out,"\n\n================ \n= HGT #%d %s \n==============",i,(diffHGT==1)?" (Trivial transfer: RF distance decreased by 1)":"");
			fprintf(res,"1");
            fprintf(res,"\nHGT %d %s",i,(diffHGT==1)?" Trivial":"");
			nbTrivial += tabHGT[i].trivial; //(diffHGT==1)?1:0;
	//		if(tabBoot != NULL)
	//			fprintf(out,"\nBootstrap value : %3.1lf%%",(double)tabBoot[i]);
			fprintf(res,"\n");
			for(int m=1;m<=tabHGT[i].listSource[0];m++){
				if(m == tabHGT[i].listSource[0])
					fprintf(res,"%s",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
				else
					fprintf(res,"%s, ",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
			}
			//fprintf(res,"\n%d",tabHGT[i].listDestination[0]);
			fprintf(res,"\n");
			for(int m=1;m<=tabHGT[i].listDestination[0];m++){
				if(m == tabHGT[i].listDestination[0])
					fprintf(res,"%s",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);
				else
					fprintf(res,"%s, ",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);
			}
	//		printTransfer(out,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
			printTransfer(res,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
	//		fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
			fprintf(res,"\nRF = %d , LS = %1.3lf , BD = %1.1lf\n",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
		}
		if(nbHGT == 0){
	//		fprintf(out,"\n\nNo more HGT have been detected due to evolutionary constraints");
	//		fprintf(res,"\nNo more HGT have been detected due to evolutionary constraints\n");
		}
		else
			if(tabHGT[nbHGT].crit.RF != 0){
	//			fprintf(out,"\n\nNo more HGT have been detected due to evolutionary constraints ");
	//			fprintf(res,"\nNo more HGT have been detected due to evolutionary constraints\n");
			}
    
	} */
	// else{
		
		i=1;
		int lastRF=0;
		int trivial;
		int k1=1;
		for(j=1;j<=multicheckTab[0].m;j++){
		
			if(multicheckTab[j].nbHgtFound == 0) continue;
			//	fprintf(out,"\n\n=======> %d %s been found for this iteration",multicheckTab[j].nbHgtFound,(multicheckTab[j].nbHgtFound==1)?"hgt has":"hgts have");
				fprintf(res,"%d",multicheckTab[j].nbHgtFound);
				if(((j == multicheckTab[0].m) && (strcmp(stepbystep,"yes") == 0)) || (strcmp(stepbystep,"yes") != 0) ){
					fprintf(res2,"%d",multicheckTab[j].nbHgtFound);
				}
			for(k=1;k<=multicheckTab[j].nbHgtFound;k++){
				if(j==1)
					diffHGT = RFref - tabHGT[i].crit.RF;
				else
					diffHGT = multicheckTab[j-1].RF - tabHGT[i].crit.RF;
				trivial=0;
				if(
				  (tabHGT[i].source_A == tabHGT[i].dest_A ) ||
				  (tabHGT[i].source_A == tabHGT[i].dest_B ) ||
				  (tabHGT[i].source_B == tabHGT[i].dest_A ) ||
				  (tabHGT[i].source_B == tabHGT[i].dest_B ) 
				){
					trivial=1;
				}
			//			fprintf(out,"\n================ \n= HGT #%d/%d %s\n================",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 1 )?" (Trivial transfer: RF distance decreased by 1)":"");//((trivial==1) && (diffHGT==1) )?" (Trivial transfer: RF distance decreased by 1)":"");
				fprintf(res,"\nHGT %d / %d %s%s%s%s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 1 )?" Trivial":"",((tabHGT[i].valide == 2 ) || (tabHGT[i].valide == 6))?" Reverse":"",(tabHGT[i].valide == 3 )?" Canceled":"",(tabHGT[i].valide == 5 )?" Initial":"");
				// fprintf(res,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].valide == 2 )?" Reverse":"");
				// fprintf(res,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].valide == 3 )?" Canceled":"");
				// fprintf(res,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].valide == 5 )?" Initial":"");
				// fprintf(res,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 6 )?" Reverse":"");
				
				if(((j == multicheckTab[0].m) && (strcmp(stepbystep,"yes") == 0)) || (strcmp(stepbystep,"yes") != 0) ){
					fprintf(res2,"\nHGT %d / %d %s%s%s%s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 1 )?" Trivial":"",((tabHGT[i].valide == 2 ) || (tabHGT[i].valide == 6))?" Reverse":"",(tabHGT[i].valide == 3 )?" Canceled":"",(tabHGT[i].valide == 5 )?" Initial":"");
					// fprintf(res2,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].valide == 2 )?" Reverse":"");
					// fprintf(res2,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].valide == 3 )?" Canceled":"");
					// fprintf(res2,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].valide == 5 )?" Initial":"");
					// fprintf(res2,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 6 )?" Reverse":"");
					fprintf(res2,"\nHGT %d / %d %s%s%s%s",k1++,nbHGT,(tabHGT[i].trivial == 1 )?" Trivial":"", ((tabHGT[i].valide == 2 ) || (tabHGT[i].valide == 6))?" Reverse":"",(tabHGT[i].valide == 3 )?" Canceled":"",(tabHGT[i].valide == 5 )?" Initial":"");
					
				}
				if(tabHGT[i].trivial == 1)
					nbTrivial ++;
			  //fprintf(out,"\n%s(%d)",(tabHGT[i].trivial == 1)?"TOTO":"TATA",tabHGT[i].trivial);
		//      if(tabBoot != NULL)
			//				fprintf(out,"\nBootstrap value : %3.1lf%%",(double)tabBoot[i]);

				if((tabBoot != NULL && tabBoot[i] >= bootmin) || (tabBoot == NULL)){
					//fprintf(res,"\n%d",tabHGT[i].listSource[0]);
					fprintf(res,"\n");
					for(int m=1;m<=tabHGT[i].listSource[0];m++){
						//fprintf(res,"%s ",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
						if(m == tabHGT[i].listSource[0])
							fprintf(res,"%s",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
						else
							fprintf(res,"%s, ",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
					}
					//fprintf(res,"\n%d",tabHGT[i].listDestination[0]);
					fprintf(res,"\n");
					for(int m=1;m<=tabHGT[i].listDestination[0];m++){
						//	fprintf(res,"%s ",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);
						if(m == tabHGT[i].listDestination[0])
							fprintf(res,"%s",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);	
						else
							fprintf(res,"%s, ",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);
								
					}
			//			printTransfer(out,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
					printTransfer(res1,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
					if(((j == multicheckTab[0].m) && (strcmp(stepbystep,"yes") == 0)) || (strcmp(stepbystep,"yes") != 0) ){
						//fprintf(res2,"\n%d--%d %d--%d",tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
						printTransfer(res2,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
					}
			//		  fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
//					fprintf(res1,"\nRF = %d , LS = %1.3lf , BD = %1.1lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
//					fprintf(res1,"\nrRF = %d , rLS = %1.3lf , rBD = %1.1lf",tabHGT[i].crit.rRF,tabHGT[i].crit.rLS,tabHGT[i].crit.rBD,tabHGT[i].crit.QD);
					if(((j == multicheckTab[0].m) && (strcmp(stepbystep,"yes") == 0)) || (strcmp(stepbystep,"yes") != 0) ){
//						fprintf(res2,"\nDIRECT : RF = %d , LS = %1.3lf , BD = %1.1lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
//						fprintf(res2,"\nREVERSE: RF = %d , LS = %1.3lf , BD = %1.1lf",tabHGT[i].crit.rRF,tabHGT[i].crit.rLS,tabHGT[i].crit.rBD,tabHGT[i].crit.QD);
					}
					i++;
				//	}
				}
			}
			//	fprintf(out,"\n\n=======> After this iteration the criteria are :");
		//		fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",multicheckTab[j].RF,multicheckTab[j].LS,multicheckTab[j].BD,multicheckTab[j].QD);
//			fprintf(res1,"\nRF = %d , LS = %1.3lf , BD = %1.1lf\n",multicheckTab[j].RF,multicheckTab[j].LS,multicheckTab[j].BD,multicheckTab[j].QD);
			if(((j == multicheckTab[0].m) && (strcmp(stepbystep,"yes") == 0)) || (strcmp(stepbystep,"yes") != 0) ){
//				fprintf(res2,"\nRF = %d , LS = %1.3lf , BD = %1.1lf\n",multicheckTab[j].RF,multicheckTab[j].LS,multicheckTab[j].BD,multicheckTab[j].QD);
			}
			lastRF = multicheckTab[j].RF;
		}
	//		fprintf(out,"\n\nTotal number of HGTs : %d",nbHGT);
			
		if(lastRF != 0){
			//	fprintf(out,"\n\nNo more HGT have been detected due to evolutionary constraints");
		//		fprintf(res,"\nNo more HGT have been detected due to evolutionary constraints\n");
		}
		
	// }
	fprintf(res1,"Stat= %d %d",nbHGT,nbTrivial);
}
//=================================================================
//==
//=================================================================
void printBestHGT_F(FILE *out,int noHGT,struct InputTree SpeciesTree,struct HGT bestHGT,int *tmp,int nbTree, int boot){
	
	if((bestHGT.valide == 1) || (bestHGT.valide == 6))
	{

		/*printf("\n\n======== Transfer %d ==========",noHGT-(*tmp));
		printf("\n%d--%d -> %d--%d",bestHGT.source_A,bestHGT.source_B,bestHGT.dest_A,bestHGT.dest_B);
		printf("\nRF = %d \nLS = %lf \nBD = %lf",bestHGT.crit.RF,bestHGT.crit.LS,bestHGT.crit.BD);*/

		fprintf(out,"\n\n================ \n= HGT #%d \n==============",noHGT-(*tmp));
		if(boot != -1)
			fprintf(out,"\n Boottrap value : %lf%%",(double)boot/nbTree*100.0);
		printTransfer(out,1,SpeciesTree.SpeciesName, SpeciesTree.size,bestHGT.source_A,bestHGT.source_B,bestHGT.dest_A,bestHGT.dest_B);
		//printf("\n%d--%d -> %d--%d",SpeciesTree.ARETE[2*(bestHGT.source)-1],SpeciesTree.ARETE[2*(bestHGT.source)-2],SpeciesTree.ARETE[2*(bestHGT.destination)-1],SpeciesTree.ARETE[2*(bestHGT.destination)-2]);
		fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",bestHGT.crit.RF,bestHGT.crit.LS,bestHGT.crit.BD,bestHGT.crit.QD);
/*
		fprintf(out,"\n%d",bestHGT.listSource[0]);
		for(j=1;j<=bestHGT.listSource[0];j++)
			fprintf(out," %d",bestHGT.listSource[j]);
		fprintf(out,"\n%d",bestHGT.listDestination[0]);
		for(j=1;j<=bestHGT.listDestination[0];j++)
			fprintf(out," %d",bestHGT.listDestination[j]);
*/	}
	else
		(*tmp) =+ 1;
}




int compareLeaves(int *listLeaves,int *tab,int n){
	int i,/*val,*/cpt1=0,cpt2=0;
	
	for(i=1;i<=n;i++){
		if(listLeaves[i] != tab[i]) cpt1++;
		tab[i] = (tab[i]+1) % 2;
	}
	for(i=1;i<=n;i++){
		if(listLeaves[i] != tab[i]) cpt2++;
	}

	if (cpt1 < cpt2) return cpt1;
	else return cpt2;

}
//============================================================================
//== 
//============================================================================
int findApproxRootBranch(int * listLeaves,struct InputTree *aTree, int maxDiff){
	
	int choix =-1,i,j,/*k,*/val,valmax = maxDiff;

	int * tab = (int *) malloc((aTree->size+1)*sizeof(int));

	//= recherche de la branche equivalente a listLeaves

	for(i=1;i<=2*aTree->size-3-aTree->kt;i++){
		for(j=1;j<=aTree->size;j++)
			if(aTree->ADD[j][aTree->ARETE[2*i-1]] < aTree->ADD[j][aTree->ARETE[2*i-2]])
				tab[j] = 0;
			else
				tab[j] = 1;
		/*printf("\n (%d)\t: ",i);
		for(int k=1;k<=aTree->size;k++){
				printf("%d ",tab[k]);
		}*/
		val = compareLeaves(listLeaves,tab,aTree->size);
		if(val <= valmax){
			valmax = val;
			choix = i;
		}
		for(int k=1;k<=aTree->size;k++){
				tab[k] = -1;
		}
	}
	//printf("\nHGT-DETECTION : no branch selectionne = %d , similarite = %d",choix,valmax);
	free(tab);

	return choix;
}

//==========================================================================
//== selectionne la branche voisine du midpoint qui minimise la distance de 
//== Robinson and Foulds 
//==========================================================================
int	bestRFNeighbor(struct InputTree *geneTree,struct InputTree *speciesTree,int position){
	int choix=position;
	int v1=-1,v2=-1,v3=-1,v4=-1,i;
	struct InputTree tmpTree;	//= arbre temporaire qui va etre modifier pour calculer RF
	struct CRITERIA aCrit;		//== struture of all the criteria
	int min= (int) INFINI;

	initInputTree(&tmpTree);	

	//== recherche des voisins
	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){

		if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] || 
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v1==-1)){
		
			v1 = i; //printf("\n(v1=%d,",i);
		}
		else if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] || 
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v2==-1)){
		
			v2 = i; //printf("v2=%d,",i);
		}
		else if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] || 
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v3==-1)){
		
			v3 = i; //printf("v2=%d,",i);
		}
		else if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] || 
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v4==-1)){
		
			v4 = i; //printf("v3=%d)",i);
		}
	}
	
	//printf("\n Selection de la racine :");

	copyInputTree(&tmpTree,(*geneTree),1,1);
	addLeafAndUpdate(&tmpTree,position);
	tmpTree.Root = tmpTree.size;
		
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

	//printf("\n(%d), %d--%d = %d",position,geneTree->ARETE[2*position-1],geneTree->ARETE[2*position-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = position;
	}
		
	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v1);
	tmpTree.Root = tmpTree.size;
		
	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);

	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

	//printf("\n(%d), %d--%d = %d",v1,geneTree->ARETE[2*v1-1],geneTree->ARETE[2*v1-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v1;
	}
	
	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v2);
	tmpTree.Root = tmpTree.size;
		
	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

	//printf("\n(%d), %d--%d = %d",v2,geneTree->ARETE[2*v2-1],geneTree->ARETE[2*v2-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v2;
	}
		
	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v3);
	tmpTree.Root = tmpTree.size;
		
	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);
	//printf("\n(%d), %d--%d = %d",v3,geneTree->ARETE[2*v3-1],geneTree->ARETE[2*v3-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v3;
	}
	
	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v4);
	tmpTree.Root = tmpTree.size;
		
	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);
	//printf("\n(%d), %d--%d = %d",v4,geneTree->ARETE[2*v4-1],geneTree->ARETE[2*v4-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v4;
	}
	
	
	FreeCriteria(&aCrit,tmpTree.size);
	return choix;
}

int compareBranches(int * tab1,int * tab2,int size){
	int score1=0,score2=0,i;
	
	for(int i=1;i<=size;i++){
		score1 = (tab1[i] != tab2[i])?1:0;
		score2 = (tab1[i] == tab2[i])?1:0;
	}
	
	return (score1<score2)?score1:score2;
}

//==========================================================================
//== selectionne la branche qui minimise la distance de Robinson and Foulds 
//==========================================================================
/*int bestbipartition(struct InputTree *geneTree, struct InputTree *speciesTree){

	int indiceBranche,i,bestScore=(int)INFINI,currentScore,bestPos=-1,j;
	struct InputTree tmpTree;		//= arbre temporaire qui va etre modifier pour calculer RF
	struct CRITERIA aCrit;			//== struture of all the criteria
	int choix=-1;
	int min= (int) INFINI;

	int * PLACEk1=(int *) malloc((2*geneTree->size-2)*sizeof(int));
	int * PLACEk2=(int *) malloc((2*geneTree->size-2)*sizeof(int));
	
	int ** Bk1=(int **) malloc((2*geneTree->size-2)*sizeof(int*));
	int ** Bk2=(int **) malloc((2*geneTree->size-2)*sizeof(int*));

	for (i=0;i<2*geneTree->size-2;i++)
	{
		Bk1[i]=(int *) malloc((geneTree->size)*sizeof(int));
		Bk2[i]=(int *) malloc((geneTree->size)*sizeof(int));     
	}
	
	Bipartition_Table(speciesTree->ADD,Bk1,PLACEk1,geneTree->size+1);
	Bipartition_Table(geneTree->ADD.Matrice,Bk2,PLACEk2,geneTree->size+1);
					
	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){
		for(j=1;j<=2*geneTree->size-3-geneTree->kt;j++){			
			currentScore = compareBranches(Bk1[i],Bk2[j],geneTree->size);
			if(currentScore < bestScore){
				bestScore = currentScore;
				bestPos = i;
				if(bestScore == 0){
					i=j=(int)INFINI;
				}
			}
		}
	}
	
	//============
	int pos_a=-1,pos_b=-1;
	int nbUn=0;
	int nbZero=0;
	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){
		if()
	}
	
	
	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){
		if(geneTree->ADD[geneTree->ARETE[2*i-1]][geneTree->ARETE[2*i-1]]){
		}
	}
	

	return choix;
}*/

//==========================================================================
//== selectionne la branche qui minimise la distance de Robinson and Foulds 
//==========================================================================
int bestRFBranch(struct InputTree *geneTree, struct InputTree *speciesTree){

	int indiceBranche;
	struct InputTree tmpTree;	//= arbre temporaire qui va etre modifier pour calculer RF
	struct CRITERIA aCrit;		//== struture of all the criteria
	int choix=-1;
	int min= (int) INFINI;

	initInputTree(&tmpTree);	

	printf("\nTaille = %d",geneTree->size);
	for(indiceBranche=1;indiceBranche<=2*geneTree->size-3-geneTree->kt;indiceBranche++){
		
		copyInputTree(&tmpTree,(*geneTree),1,1);
		printf("%d",indiceBranche);
		addLeafAndUpdate(&tmpTree,indiceBranche);
		exit(0);
		tmpTree.Root = tmpTree.size;
		
		InitCriteria(&aCrit,tmpTree.size);
		computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

		printf("\n(%d), %ld--%ld = %d",indiceBranche,geneTree->ARETE[2*indiceBranche-1],geneTree->ARETE[2*indiceBranche-2],aCrit.RF);
		if(aCrit.RF < min){
			min=aCrit.RF;
			choix = indiceBranche;
		}
		FreeCriteria(&aCrit,tmpTree.size);
	}

	return choix;
}

bool file_exists(const char * filename)
{
	FILE *file;

    if ((file=fopen(filename, "r"))==NULL)
    {
		return false;
    }
	else{
        fclose(file);
        return true;
	}
}

//==========================================================================================================
//== 
//==========================================================================================================
int addRoot(struct InputTree *aTree,struct InputTree *refTree,const char * message, const char * add,char *fichier,char *fichier_leaves,int * listLeaves,const char * version){
	
	int choix=-1,choix2=-1;
	FILE *in;
	int R1,R2,i,j;
	
	if(strcmp(add,"prompt") == 0 ){
		printf("\n%s",message);
		printEdges(NULL,1,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,aTree->size,NULL,0,aTree->kt);
		
		while(choix<1 || choix > 2*aTree->size-3-aTree->kt){
			printf("\n\nSelect the root branch :");
			scanf("%d",&choix2);
		}
	}
	else if(strcmp(add,"file") == 0 ){
		
		if(file_exists(fichier_leaves)){
			int * listLeaves2 = (int*) malloc((aTree->size+1)*sizeof(int));
			/*if(listLeaves2 == NULL){
				printf("\nlistLeaves2 == NULL");
			}*/
			for(i=1;i<=aTree->size;i++)
				listLeaves2[i] = 2;
			
			if((in=fopen(fichier_leaves,"r"))==NULL){
				printf("\nHGT-DETECTION : Can't open root file (%s)",fichier_leaves);
				exit(-1);
			}
	
			char * tmp;
			tmp = (char*)malloc(100);
			//printf("\nHGT-DETECTION (size = %d): ",aTree->size);
			
			do{
				fscanf(in,"%s",tmp);
				//printf("\n%s ",tmp);
				if(strcmp(tmp,"<>") != 0){
					for(i=1;i<=aTree->size;i++){
						if(strcmp(aTree->SpeciesName[i],tmp) == 0){
							listLeaves2[i] = 0;
							//printf(" (0) ");
						}
					}
				}			
			}while(strcmp(tmp,"<>") != 0);
			
			while(fscanf(in,"%s",tmp) != -1){
				//printf("\n%s ",tmp);
				if(strcmp(tmp,"<>") != 0){
					for(i=1;i<=aTree->size;i++){
						if(strcmp(aTree->SpeciesName[i],tmp) == 0){
							listLeaves2[i] = 1;
							//printf(" (1) ");
						}
					}	
				}
			}
			
			fclose(in);
			
			/*printf("\n");
			for(i=1;i<=aTree->size;i++){
				printf(">%s<",aTree->SpeciesName[i]);
			}	
			printf("\nHGT-DETECTION : \n\t\t  ");
			for(i=1;i<=aTree->size;i++){
				printf("%d ",listLeaves2[i]);
			}	*/
			
			choix = choix2 = findApproxRootBranch(listLeaves2,aTree,10000);
			
			// fscanf(in,"%d %d",&R1,&R2);
			// for(i=1;i<=2*aTree->size-3-aTree->kt;i++){
				// if((aTree->ARETE[2*i-1] == R1 && aTree->ARETE[2*i-2] == R2) || (aTree->ARETE[2*i-1] == R2 && aTree->ARETE[2*i-2] == R1))
					// choix2 = choix = i;
			// }
			
			//printf("\n\n====%d====\n",choix);
			
		
		}else{
			if((in=fopen(fichier,"r"))==NULL){
				printf("\nCan't open root file (%s)",fichier);
				exit(-1);
			}
			fscanf(in,"%d %d",&R1,&R2);
			for(i=1;i<=2*aTree->size-3-aTree->kt;i++){
				if((aTree->ARETE[2*i-1] == R1 && aTree->ARETE[2*i-2] == R2) || (aTree->ARETE[2*i-1] == R2 && aTree->ARETE[2*i-2] == R1))
					choix2 = choix = i;
			}
			fclose(in);
		}
	}
	else if(strcmp(add,"bestrfbranch") == 0){
		printf("choix = %d",choix);	
		choix2 = choix = bestRFBranch(aTree,refTree);
		printf("choix = %d",choix);	
		//listLeaves = NULL;
	}
	/*else if(strcmp(add,"bestbipartition") == 0){
		printf("choix = %d",choix);	
		choix2 = choix = bestbipartition(aTree,refTree);
		printf("choix = %d",choix);	
		//listLeaves = NULL;
	}*/
	else{
		choix2 = choix = midPoint(aTree->ARETE,aTree->ADD,aTree->size,aTree->kt);
		if(refTree != NULL){
			choix2 = choix = bestRFBranch(aTree,refTree);
			//choix2 = choix = bestRFNeighbor(aTree,refTree,choix);
		}
	}
	
	//printf("\nchoix = %d",choix);
	if(listLeaves != NULL){
		//== sauvegarder la liste des feuilles de chaque cot� de l'arete.
		if(listLeaves[0] == -1){
			for(i=1;i<=aTree->size;i++)
				if(aTree->ADD[i][aTree->ARETE[2*choix-1]] < aTree->ADD[i][aTree->ARETE[2*choix-2]])
					listLeaves[i] = 0;
				else
					listLeaves[i] = 1;
			listLeaves[0]=1;
		}
		else{
			//printf("\nHGT-DETECTION : Recherche de la racine en fonction de la bipartition");
			if(strcmp(add,"bestbipartition") == 0){
				choix2 = findApproxRootBranch(listLeaves,aTree,10000);
			}
			//printf("\nHGT-DETECTION : choix = %d",choix2);
		}

		if(choix2==-1)
			return -1;
		else
			choix=choix2;
	}

	printRoot(fichier,aTree->ARETE[2*choix-1],aTree->ARETE[2*choix-2]);
	printRootByLeaves(fichier_leaves,choix,aTree);

	addLeafAndUpdate(aTree,choix);
	aTree->Root = aTree->size;
	
	//===mise a jour de degre
	for(i=1;i<=2*(aTree->size)-2-aTree->kt;i++){
		aTree->degre[i]=0;
		for(j=1;j<=2*aTree->size-2-aTree->kt;j++)
			if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
	}
	return 0;

}


//=================================================================
//== 
//=================================================================
void chargerFichier(struct InputTree *aTree,const char* branchesFile,const char* rootFile){
	
	FILE *branch = fopen(branchesFile,"r");
	FILE *root =  fopen(rootFile,"r");
	int aa,bb,i,noBranch=0;
	double ll;

	for(i=1;i<=2*(aTree->size)-3;i++){
		fscanf(branch,"%d %d %lf",&aa,&bb,&ll);
		aTree->ARETE[2*i-1] = aa;
		aTree->ARETE[2*i-2] = bb;
		aTree->LONGUEUR[i-1] = ll;
		//printf("\n%d %d %lf",ARETE[2*i-1],ARETE[2*i-2],LONGUEUR[i-1]);
	}
	fscanf(root,"%d %d",&aa,&bb);

	for(i=1;i<=2*(aTree->size)-3;i++)
		if((aTree->ARETE[2*i-1] == aa && aTree->ARETE[2*i-2] == bb) || (aTree->ARETE[2*i-1] == bb && aTree->ARETE[2*i-2] == aa))
			noBranch = i;
	printf("-->nobranch=%d<--",i);
	addLeafAndUpdate(aTree,noBranch);
	aTree->Root = aTree->size;

	fclose(branch);
	fclose(root);
}


//=================================================================
//== read the species or the gene tree in the input file
//=================================================================
int readInput(int Type, const char *file,struct InputTree * aTree){

	int size,i,j;
	char name[50];
	double val;
	FILE * in;
	
	//= ouverture du fichier
	if((in = fopen(file,"r"))== NULL)
		return -1;

	//= lecture de la taille des matrices
	fscanf(in,"%d",&size);
	
	//= allocation de la m�moire
	//allocMemmory(aTree,1);
	aTree->size = size;
	size++; // more space for the root
	aTree->SpeciesName = (char **)malloc((size+1)*sizeof(char*));
	aTree->Input = (double**)malloc((2*size)*sizeof(double*));
	aTree->ADD = (double**)malloc((2*size)*sizeof(double*));
	aTree->W = (double**)malloc((size+1)*sizeof(double*));
	for(i=0;i<2*size;i++){
		aTree->ADD[i] = (double*)malloc((2*size)*sizeof(double));
		aTree->Input[i] = (double*)malloc((2*size)*sizeof(double));
		if(i<=size){
			aTree->SpeciesName[i] = (char*)malloc(SPECIES_NAME_LENGTH);
			aTree->W[i] = (double*)malloc((size+1)*sizeof(double));
		}
	}

	for(i=0;i<=size;i++)
		for(j=0;j<=size;j++)
			aTree->W[i][j] = 1.0;

	size--;
	//= reads species tree
	for(i=1;i<=size;i++)
	{
		fscanf(in,"%s",name);
		if(Type == SPECIE) strcpy(aTree->SpeciesName[i],name);
		for(j=1;j<=size;j++)
		{
			fscanf(in,"%lf",&val);
			if(Type == SPECIE) aTree->Input[i][j] = val;
		}	    
	}

	//= read gene tree
	for(i=1;i<=size;i++)
	{
		fscanf(in,"%s",name);
		if(Type == GENE) strcpy(aTree->SpeciesName[i],name);
		for(j=1;j<=size;j++)
		{
			fscanf(in,"%lf",&val);
			if(Type == GENE) aTree->Input[i][j] = val;
		}	
	}

	strcpy(aTree->SpeciesName[size+1],"Root");
	
	fclose(in);
	
	return 0;
}
//===========================================================
//== print the input matrix
//===========================================================
void printMatrix(char** Name, double **Matrix,int size){

	int i,j;
	
	for(i=1;i<=size;i++){
//		printf("\n%s",Name[i]);
		printf("\n%d",i);
		for(j=1;j<=size;j++){
			if(Matrix[i][j] >= INFINI)
				printf("  INFINI ");
			else
				printf(" %3.5lf",Matrix[i][j]);
		}
	}
}
//===========================================================
//== print the input matrix
//===========================================================
void printMatrix2(char** Name, double **Matrix,int size){

	int i,j;
	
	for(i=1;i<=size;i++){
		//printf("\n%s",Name[i]);
		printf("\n%d",i);
		for(j=1;j<=size;j++){
			if(Matrix[i][j] < epsilon )
				printf("  NULL ");
			else
				printf(" %3.5lf",Matrix[i][j]);
		}
	}
}

//==============================================================
//== print branches
//==============================================================
void printBranches(FILE *out,struct InputTree aTree,const char * Message,int * BSARETE,int nbTree){

	if(out==NULL)
		printf("%s",Message);
	else
		fprintf(out,"\n\n%s",Message);
	printEdges(out,1,aTree.ARETE,aTree.LONGUEUR,aTree.SpeciesName,aTree.size,BSARETE,nbTree,aTree.kt);
}



//=====================================================================
//
//=====================================================================
int ExtraireDonnees(const char * chaine, char *champs, char * contenu){

	int cpt=0,i;
	int egale=false;
	int tailleChaine;

	if(chaine[0] != '-')
		return false;
	tailleChaine = (int)strlen(chaine);

	for(i=1;i<tailleChaine;i++){

		if (chaine[i] == '='){
			egale = true;
			champs[cpt] = '\0';
			cpt=0;
			continue;
		}
		if(egale)
			contenu[cpt++] = chaine[i];
		else
			champs[cpt++] = chaine[i];
	}
	contenu[cpt] = '\0';

	if (!egale)
		return false;

	return true;
}
//===================================================================================
//=
//===================================================================================
int readParameters(struct Parameters * param){

	char champs[100];
	char contenu[100];
	char input[100];
	char output[100];
	char hgtResultFile[100];
	int i;
	sprintf((*param).sort,"yes");
  	sprintf((*param).printWeb,"yes");
	sprintf((*param).criterion,"bd");
	sprintf((*param).verbose,"no");
	sprintf((*param).mode,"multicheck");
	sprintf((*param).viewtree,"no");
	sprintf((*param).generoot,"midpoint");
	sprintf((*param).speciesroot,"midpoint");
	sprintf((*param).load,"no");
	sprintf((*param).version,"consol");
	sprintf((*param).multiple,"no");
	sprintf((*param).multigene,"no");
	sprintf((*param).path,"");
	sprintf(input,"");
	sprintf(output,"output.txt");
	sprintf(hgtResultFile,"hgtresultfile.txt");
	sprintf((*param).scenario,"unique");
	sprintf((*param).subtree,"yes");
	sprintf((*param).bootstrap,"no");
	(*param).nbhgt = 100;
	(*param).bootmin = 0;
	sprintf((*param).special,"no");
	(*param).rand_bootstrap = 0;
	sprintf((*param).stepbystep,"no");
	

	sprintf((*param).inputfile,"%s%s",(*param).path,input);
	sprintf((*param).input,"%sinput_.txt",(*param).path);
	sprintf((*param).outputfile,"%s%s",(*param).path,output);
	sprintf((*param).results,"%sresults.txt",(*param).path);
	sprintf((*param).results_bouba,"%sresults2.txt",(*param).path);
	sprintf((*param).hgtResultFile,"%s%s",(*param).path,hgtResultFile);
	sprintf((*param).speciesTree,"%sspeciesTree.txt",(*param).path);
	sprintf((*param).geneTree,"%sgeneTree.txt",(*param).path);
	sprintf((*param).speciesRootfile,"%sspeciesRoot.txt",(*param).path);
	sprintf((*param).speciesRootfileLeaves,"%sspeciesRootLeaves.txt",(*param).path);
	sprintf((*param).geneRootfile,"%sgeneRoot.txt",(*param).path);
	sprintf((*param).errorFile,"%serrorFile.txt",(*param).path);
	sprintf((*param).geneRootfileLeaves,"%sgeneRootLeaves.txt",(*param).path);
	sprintf((*param).speciesTreeWeb,"%sspeciesTreeWeb.txt",(*param).path);
	sprintf((*param).geneTreeWeb,"%sgeneTreeWeb.txt",(*param).path);
	sprintf((*param).outputWeb,"%soutputWeb.txt",(*param).path);
	sprintf((*param).noMoreHgtfile,"%snomorehgt.txt",(*param).path);
	sprintf((*param).prehgtfile,"%sprehgt.txt",(*param).path);
	
	return 0;
}

void saveTree(char * fichier,struct InputTree SpeciesTree,struct HGT *bestHGT, int addHGT,int cpt_hgt,const char *subtree,char *SCENARIO,double *bootStrap){
	
	FILE *out;
	int i, first=0;

	if ((out=fopen(fichier,"w+"))==NULL){
		printf("can't open species web file");
		exit(-1);
	}

	fprintf(out,"maptype = VERTICAL");
	fprintf(out,"\nintVerPrint = O");
	fprintf(out,"\nproportion = N");
	fprintf(out,"\ndrawRealNames = 1");
	fprintf(out,"\ncolorobject = BLACK");
	fprintf(out,"\ncolorreticulation = MAGENTA");
	fprintf(out,"\ncoloredge = BLUE");
	fprintf(out,"\nkt = %d",SpeciesTree.kt);
	fprintf(out,"\ndrawRetic = %s",(addHGT == 1)?"2":"0");
	fprintf(out,"\nn = %d",SpeciesTree.size);
	fprintf(out,"\net = ");
	for(i=1;i<=SpeciesTree.size;i++){
		if(i!=1) fprintf(out,",");
		fprintf(out,"%s",SpeciesTree.SpeciesName[i]);
	}
	fprintf(out,"\naretes = ");
	for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
		if(i!=1) fprintf(out,",");
		fprintf(out,"%ld,%ld",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2]);
	}
	fprintf(out,"\nlongueur = ");
	for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
		if(i!=1) fprintf(out,",");
		fprintf(out,"%lf",SpeciesTree.LONGUEUR[i-1]);
	}
	if(addHGT == 1){
		fprintf(out,"\nhgt = ");
		for(i=1;i<=cpt_hgt;i++){
			if((bestHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/){
				if(first==1) fprintf(out,","); 
				first=1;
				fprintf(out,"%d,%d,%d,%d",bestHGT[i].source_A,bestHGT[i].source_B,bestHGT[i].dest_A,bestHGT[i].dest_B);
			}
		}
		fprintf(out,"\nhgt_reels = ");
		first=0;
		for(i=1;i<=cpt_hgt;i++){
			if((bestHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/){
				if(first==1) fprintf(out,","); 
				first=1;
				fprintf(out,"%d,%d,%d,%d",bestHGT[i].source_Ar,bestHGT[i].source_Br,bestHGT[i].dest_Ar,bestHGT[i].dest_Br);
			}
		}
		if(bootStrap != NULL){
			first=0;
			fprintf(out,"\nbootHGT = ");
			for(i=1;i<=cpt_hgt;i++){
				if(first==1) fprintf(out,","); 
				first=1;
				fprintf(out,"%1.0lf",bootStrap[i]);
			}
		}
		fprintf(out,"\nhgt_valide = ");
		first=0;
		for(i=1;i<=cpt_hgt;i++){
			if((bestHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/){
				if(first==1) fprintf(out,","); 
				first=1;
				
				fprintf(out,"%d",bestHGT[i].valide);
			}
		}
		fprintf(out,"\nhgt_sequence = ");
		first=0;
		for(i=1;i<=cpt_hgt;i++){
			if((bestHGT[i].valide >=1)/* && (bestHGT[i].trivial==0)*/){
				if(first==1) fprintf(out,","); 
				first=1;
				fprintf(out,"%d",bestHGT[i].sequence);
			}
		}
	}
	if(addHGT == 1)
		fprintf(out,"\nroot=%d",SpeciesTree.Root);
	else
		fprintf(out,"\nroot=0");

	fclose(out);
}



//=============================================================
//
//=============================================================
void deleteSommet(struct DescTree *tab, int t1, int t2){

	int i,j,k=0,toDel=1,tmp;
	
	// si les elements de tab1 se trouve dans tab2 on les supprimes dans tab2
	for(i=0;i<tab[t1].nbSommet;i++){
		tmp=0;
		for(j=0;j<tab[t2].nbSommet;j++){
			if(tab[t1].Tableau[i] == tab[t2].Tableau[j]){
				//printf("\n=>%d (%d)",tab[t1].Tableau[i],tab[t1].nbSommet);
				tmp=1;
			}
		}
		if(tmp==0) toDel=0;
	}
	for(i=0;i<tab[t1-1].nbSommet;i++){
		tmp=0;
		for(j=0;j<tab[t2].nbSommet;j++){
			if(tab[t1-1].Tableau[i] == tab[t2].Tableau[j]){
				//printf("\n=>%d (%d)",tab[t1-1].Tableau[i],tab[t1-1].nbSommet);
				tmp=1;
			}
		}
		if(tmp==0) toDel=0;
	}
	

	// suppression des sommets
	if(toDel == 1){
		for(i=0;i<tab[t1].nbSommet;i++){
			for(j=0;j<tab[t2].nbSommet;j++){
				if(tab[t1].Tableau[i] == tab[t2].Tableau[j]){
					tab[t2].Tableau[j]=-1;
					k=j;
				}
			}

			for(j=k;j<tab[t2].nbSommet-1;j++){
				tab[t2].Tableau[j] = tab[t2].Tableau[j+1];
			}
			tab[t2].nbSommet--;
		
		}
	}
}


//===================================================================
//
//===================================================================
void supprimerSousEnsemble(struct DescTree *tabDetect,int nbTableau){
	
	int i,j/*,k*/;
	
	for(i=1;i<=nbTableau;i=i+2){
		for(j=i+1;j<=nbTableau;j++){
			deleteSommet(tabDetect,i,j);
		}	
	}
}

//==========================================================================================================================================
//
//==========================================================================================================================================
int formatResult(struct HGT *tabHGT,int nbHGT, struct HGT *outHGT,struct InputTree aTree){

	int nbSommet,i,j,nbTrans,sommet,tmp=1,cpt=-1,nbHGTRet=0;
	int racine;
	double max;
	int nbArete = 2*(aTree.size)-3-aTree.kt,sequence=1,sequence_old=0;

	struct DescTree *tabDetect  = (struct DescTree*)malloc((2*(nbHGT+1))*sizeof(struct DescTree));

	//res = fopen(result,"r");

	//==lecture des transferts

	//printf("\nHGT-DETECTION : nombre de HGT avant formatRes %d",nbHGT);
	
	for(i=1;i<=nbHGT;i++){
		
		if(tabHGT[i].valide >= 1){
			nbHGTRet++;
			/*sequence_old = tabHGT[i].sequence;
			
			outHGT[nbHGTRet].valide = 1;
			if(tabHGT[i].sequence != tabHGT[i-1].sequence)
					sequence ++;
			}*/
			
			if(sequence_old == 0){
				sequence = 1;
			}else{
				if(tabHGT[i].sequence != sequence_old){
					sequence ++;
				}
			}
			sequence_old = tabHGT[i].sequence;
			outHGT[nbHGTRet].sequence = sequence;
			
			copyHGT(tabHGT[i],&outHGT[nbHGTRet]);
			cpt++;
			nbSommet = tabDetect[cpt].nbSommet = tabHGT[i].listSource[0];
			
			//printf("\ni=%d",i);
			tabDetect[cpt].Tableau = (int*)malloc(nbSommet*sizeof(int));
			//printf("\n");
			for(j=0;j<nbSommet;j++){
				tabDetect[cpt].Tableau[j] = tabHGT[i].listSource[j+1];
			//	printf("%d ",tabHGT[i].listSource[j+1]);
			}
			//printf("\ni=%d",i);
			//if(i==73) {printf("\nOn sort nbSommet=%d",nbSommet); exit(1);}
			cpt++;
			nbSommet = tabDetect[cpt].nbSommet = tabHGT[i].listDestination[0];
			tabDetect[cpt].Tableau = (int*)malloc(nbSommet*sizeof(int));
			//printf("\n");
			for(j=0;j<nbSommet;j++){
				tabDetect[cpt].Tableau[j] = tabHGT[i].listDestination[j+1];
			}
		}
		else{
			printf("\nhgt : Transfer %d was found to be idle",i);
		}
	}

	outHGT[nbHGTRet+1].valide=0;
	for(int i=1;i<=nbHGTRet;i++){
		outHGT[i].source_Ar = outHGT[i].source_A;
		outHGT[i].source_Br = outHGT[i].source_B;
		outHGT[i].dest_Ar = outHGT[i].dest_A;
		outHGT[i].dest_Br = outHGT[i].dest_B;
	}
	
	nbTrans = cpt;

	supprimerSousEnsemble(tabDetect,nbTrans);

	racine = aTree.Root;
	

	for(i=0;i<=nbTrans;i++){
		max=INFINI;
		sommet=tabDetect[i].Tableau[0];
		if(tabDetect[i].nbSommet > 1){
			// calcul de la distance entre la racine et le point d'intersection de tous les sommets
			for(j=1;j<tabDetect[i].nbSommet;j++){
			
				if(max > ( (aTree.ADD[racine][sommet] + aTree.ADD[racine][tabDetect[i].Tableau[j]] - aTree.ADD[sommet][tabDetect[i].Tableau[j]]) / 2.0 ) )
					max = (aTree.ADD[racine][sommet] + aTree.ADD[racine][tabDetect[i].Tableau[j]] - aTree.ADD[sommet][tabDetect[i].Tableau[j]]) / 2.0 ;  
			}
		
			// recherche du sommet
			for(j=1;j<=nbArete;j++){
				if( ((fabs(max - aTree.ADD[racine][aTree.ARETE[2*j-1]]) < 0.00001) && (aTree.ADD[racine][aTree.ARETE[2*j-2]] < aTree.ADD[racine][aTree.ARETE[2*j-1]]) && fabs(aTree.ADD[racine][sommet] - aTree.ADD[racine][aTree.ARETE[2*j-1]] - aTree.ADD[sommet][aTree.ARETE[2*j-1]])< 0.00001) ||
					((fabs(max - aTree.ADD[racine][aTree.ARETE[2*j-2]]) < 0.00001) && (aTree.ADD[racine][aTree.ARETE[2*j-1]] < aTree.ADD[racine][aTree.ARETE[2*j-2]]) && fabs(aTree.ADD[racine][sommet] - aTree.ADD[racine][aTree.ARETE[2*j-2]] - aTree.ADD[sommet][aTree.ARETE[2*j-2]])< 0.00001)  ) {
				
					if( fabs(aTree.ADD[racine][sommet] - aTree.ADD[racine][aTree.ARETE[2*j-1]] - aTree.ADD[aTree.ARETE[2*j-1]][sommet]) < 0.00001 ){		
						if(i % 2 == 0){
							outHGT[tmp].source_A = aTree.ARETE[2*j-1];
							outHGT[tmp].source_B = aTree.ARETE[2*j-2];
						}
						else{
							outHGT[tmp].dest_A = aTree.ARETE[2*j-1];
							outHGT[tmp].dest_B = aTree.ARETE[2*j-2];
							tmp++;
						}
					}
				}
			}
		}
		else{
			for(j=1;j<=nbArete;j++){
				if((aTree.ARETE[2*j-1] == sommet) || (aTree.ARETE[2*j-2] == sommet)){
					if(i % 2 == 0){
						outHGT[tmp].source_A = aTree.ARETE[2*j-1];
						outHGT[tmp].source_B = aTree.ARETE[2*j-2];
					}
					else{
						outHGT[tmp].dest_A = aTree.ARETE[2*j-1];
						outHGT[tmp].dest_B = aTree.ARETE[2*j-2];
						tmp++;
					}
				}
					
			}
		}
	}

	for(i=0;i<=cpt;i++){
		free(tabDetect[i].Tableau);
	}
	free(tabDetect);
	
	return nbHGTRet;
}



