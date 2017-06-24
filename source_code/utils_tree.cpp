//==============================================
//=
//==============================================
double BipartitionDistance (int **B, int ** B1,int n)
{
	int i,j,k,k1,cpt1,cpt2,cpt3,cpt4,t1,t2;
	double c1,c2,min1,min2,cptB=0.0,cptB1=0.0;
	int *flag,*flag1,**Bi,**Bi1;
	Bi=(int **) malloc((2*n-2)*sizeof(int*));
	Bi1=(int **) malloc((2*n-2)*sizeof(int*));
	for (i=0;i<2*n-2;i++)
	{
		Bi[i]=(int *) malloc((2*n)*sizeof(int));
		Bi1[i]=(int *) malloc((2*n)*sizeof(int)); 
		if ((Bi1[i]==NULL)||(Bi[i]==NULL))
		{
			printf(" Data matrix is too large\n "); 
			exit(2);
		}        
	}
	flag = (int*) malloc(2*n*(sizeof(int)));
	flag1 = (int*) malloc(2*n*(sizeof(int)));
	for(i=1;i<=2*n-3;i++){
		flag[i] = flag1[i] = 0;
		for(j=1;j<=n;j++){
			flag[i] = flag[i] + B[i][j];
			flag1[i] = flag1[i] + B1[i][j];
		}
		if(flag[i] == 1 || flag[i] == n-1)
			flag[i] = 1;
		else
			flag[i] = 0;
		if(flag1[i] == 1 || flag1[i] == n-1)
			flag1[i] = 1;
		else
			flag1[i] = 0;
	}
    k=k1=1;
	for(i=1;i<=2*n-3;i++){
		if(flag[i] == 0){
			for(j=1;j<=n;j++)
				Bi[k][j] = B[i][j];
			k++;
		}
		if(flag1[i] == 0){
			for(j=1;j<=n;j++)
				Bi1[k1][j] = B1[i][j];
			k1++;
		}
	
	}
	
	/*for(i=1;i<=2*n-3;i++)
		printf("%d ",flag[i]);
	printf("\n");
	for(i=1;i<=2*n-3;i++)
		printf("%d ",flag1[i]);
	printf("\n");*/
	for(i=1;i<=n-3;i++){
		t1=t2=0;
		for(j=1;j<=n-3;j++){
			cpt1=cpt2=cpt3=cpt4=n;
			//p=p1=0;
			for(k=1;k<=n;k++){
				/*if(flag[i] == 0 && flag1[j] == 0){
					p=1;*/
					if(Bi[i][k]==Bi1[j][k])
						cpt1--;
					else
						cpt2--;
				/*}
				if(flag1[i] == 0 && flag[j] == 0){
					p1=1;*/
					if(Bi1[i][k]==Bi[j][k])
							cpt3--;
						else
							cpt4--;
				//}
			}
			/*if(p==1){*/
			c1 = ((double)cpt1<(double)cpt2)?(double)cpt1:(double)cpt2;
				if(j==1 || c1 < min1)
					min1 = c1;
			/*}
			else
				min1=0;
			if(p1==1){*/
				c2 = ((double)cpt3<(double)cpt4)?(double)cpt3:(double)cpt4;
				if(j==1 || c2 < min2)
					min2 = c2;
			/*}
			else
				min2=0;*/
		}
		cptB  += min1;
		cptB1 += min2;
	}

	for (i=0;i<2*n-2;i++)
	{
		free(Bi[i]); free(Bi1[i]);
	}

	free(Bi);
	free(Bi1);
	free(flag); free(flag1);

	return (cptB + cptB1) / 2.0;
}

/****************************************************
* application de la methode NJ pour construire
* une matrice additive
****************************************************/
void NJ(double **D1,double **DA,int n)
{
	double **D,*T1,*S,*LP,Som,Smin,Sij,L,Lii,Ljj,l1,l2,l3;
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

		if (Lii<2*epsilon) Lii=2*epsilon;
		if (Ljj<2*epsilon) Ljj=2*epsilon;
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
	if (l1<2*epsilon) l1=2*epsilon;
	if (l2<2*epsilon) l2=2*epsilon;
	if (l3<2*epsilon) l3=2*epsilon;
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

/*/===========================================
//
//===========================================*/
void TrierTableau(int tableau[],int taille){

	int i, inversion,tmp;
	do
	{
		inversion=0;

		for(i=1;i<taille;i++)
		{
			if (tableau[i]>tableau[i+1])
			{
				tmp = tableau[i];
				tableau[i] = tableau[i+1];
				tableau[i+1] = tmp;
				inversion=1;
			}
		}
	}
	while(inversion);
}

/*=======================================
//
//=======================================*/
int vecteursEgaux(struct DescTree DTSpecies, struct DescTree DTGene){

	int Tspecies = DTSpecies.nbSommet;
	int Tgene = DTGene.nbSommet;
	int i;

	if(Tspecies != Tgene) return 0;

	for(i=1;i<=Tspecies;i++)
		if(DTSpecies.Tableau[i] != DTGene.Tableau[i])
			return 0;
	return 1;
}


//===========================================================
//== Floyd
//===========================================================
void Floyd(double ** Adjacence , double ** DIST,int n,int kt)
{
	int i,j,k;

	for(i=1;i<=2*n-2-kt;i++)
		for(j=1;j<=2*n-2-kt;j++)
		{
			if(i==j)
				DIST[i][j] = 0;
			else
				DIST[i][j] = Adjacence[i][j];
		}

		for(i=1;i<=2*n-2-kt;i++)
			for(j=1;j<=2*n-2-kt;j++)
				for(k=1;k<=2*n-2-kt;k++)
				{
					if((DIST[j][i] + DIST[i][k]) < DIST[j][k])
						DIST[j][k] = DIST[j][i] + DIST[i][k];
				}
}

void Floyd(double ** Adjacence , double ** DIST,double **DIST2,int n,int kt)
{
	int i,j,k1;

	for(i=1;i<=2*n-2-kt;i++)
		for(j=1;j<=2*n-2-kt;j++)
		{
			if(i==j)
				DIST[i][j] = DIST2[i][j] = 0;
			else
				DIST[i][j] = DIST2[i][j] = Adjacence[i][j];
		}

		for(i=1;i<=2*n-2-kt;i++)
			for(j=1;j<=2*n-2-kt;j++)
				for(k1=1;k1<=2*n-2-kt;k1++)
				{
					if((DIST[j][i] + DIST[i][k1]) < DIST[j][k1])
						DIST[j][k1] = DIST2[j][k1] = DIST[j][i] + DIST[i][k1];
				}
}
//===========================================================
//== print branches
//===========================================================
void printEdges(FILE *out,int ShowName, long int * ARETES,double * LONGUEUR,char ** NAMES,int size,int *BSARETE,int nbTree,int kt){

	int i,j,fichier=1,boot=1,longChaine=0;
//	double bsvalue;

	//if(nbTree==0)
	//	boot=0;
	for(i=1;i<=size;i++){
		if ((int)strlen(NAMES[i]) > longChaine)
			longChaine = (int)strlen(NAMES[i]);
	}
	if(out == NULL)
		fichier =0;
	//printf("\n");
	for(i=1;i<=2*size-3-kt;i++){
		//if(boot==1)
			//bsvalue = BSARETE[i]*100.0 / nbTree;
		if(ShowName == 0){
			if(fichier==1){
				fprintf(out,"\n%d. %ld--%ld \t: %lf ",i,ARETES[2*i-1],ARETES[2*i-2],LONGUEUR[i-1]);
			}
			else{
				printf("\n%d. %ld--%ld \t: %lf ",i,ARETES[2*i-1],ARETES[2*i-2],LONGUEUR[i-1]);
			}
		}
		else{
			if(ARETES[2*i-1] <= size){
				if(fichier==1){
					fprintf(out,"\n%d. %ld--%s",i,ARETES[2*i-2],NAMES[ARETES[2*i-1]]);
					for(j=0;j<longChaine - (int)strlen(NAMES[ARETES[2*i-1]]);j++) fprintf(out," ");
					fprintf(out,"\t%lf",LONGUEUR[i-1]);
				}
				else{
					printf("\n%d. %ld--%s",i,ARETES[2*i-2],NAMES[ARETES[2*i-1]]);
					for(j=0;j<longChaine - (int)strlen(NAMES[ARETES[2*i-1]]);j++) printf(" ");
					printf("\t%lf",LONGUEUR[i-1]);
				}
			}
			else if(ARETES[2*i-2] <= size){
				if (fichier==1){
					fprintf(out,"\n%d. %ld--%s",i,ARETES[2*i-1],NAMES[ARETES[2*i-2]]);
					for(j=0;j<longChaine - (int)strlen(NAMES[ARETES[2*i-2]]);j++) fprintf(out," ");
					fprintf(out,"\t%lf",LONGUEUR[i-1]);
				}
				else{
					printf("\n%d. %ld--%s",i,ARETES[2*i-1],NAMES[ARETES[2*i-2]]);
					for(j=0;j<longChaine - (int)strlen(NAMES[ARETES[2*i-2]]);j++) printf(" ");
					printf("\t%lf",LONGUEUR[i-1]);
				}
			}
			else{
				if(fichier==1){
					fprintf(out,"\n%d. %ld--%ld",i,ARETES[2*i-1],ARETES[2*i-2]);
					for(j=0;j<longChaine - 2;j++) fprintf(out," ");
					fprintf(out,"\t%lf",LONGUEUR[i-1]);
				}
				else{
					printf("\n%d. %ld--%ld",i,ARETES[2*i-1],ARETES[2*i-2]);
					for(j=0;j<longChaine - 2;j++) printf(" ");
					printf("\t%lf",LONGUEUR[i-1]);
				}
			}
		}

	}
}

//=============================================================================================================
//==
//=============================================================================================================
void loadAdjacenceMatrix( double **Adjacence, long int *ARETE, double *LONGUEUR,int size,int kt){
	
	int i,j;
	
	for(i=1;i<=2*size-2;i++) /*/(n+1)*/
		for(j=1;j<=2*size-2;j++){
			Adjacence[i][j] = Adjacence[j][i] = INFINI;
      //printf("i=%d,j=%d",i,j);
}
	for(i=1;i<=2*size-3-kt;i++){
		//printf("\n%d - %d - %d",i,size,kt);
		Adjacence[ARETE[2*i-2]][ARETE[2*i-1]] = LONGUEUR[i-1];//(LONGUEUR[i-1]>5*epsilon)?LONGUEUR[i-1]:5*epsilon;
		Adjacence[ARETE[2*i-1]][ARETE[2*i-2]] = LONGUEUR[i-1]; //(LONGUEUR[i-1]>5*epsilon)?LONGUEUR[i-1]:5*epsilon;
	}
	//printf("\nhouston...");
}

//=====================================================
//==
//=====================================================
void odp1(double **D, int *X, int *i1, int *j1, int n)
{
	double S1,S;
	int i,j,k,a,*Y1;

	Y1=(int *) malloc((n+1)*sizeof(int));

	for(i=1;i<=n;i++)
		Y1[i]=1;

	X[1]=*i1;
	X[n]=*j1;
	if (n==2){
		free(Y1);
		return;
	}
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

//================================================================================
//==
//================================================================================
int Tree_edges (double **DI, long int *ARETE, double *LONGUEUR, int n,int binaire)
{

	struct EDGE { unsigned int U; unsigned int V; double LN;};
	struct EDGE *Path,*Tree;
	int i,j,k,p,P,*X;
	double S,DIS,DIS1,*L,**D;
	int pasfini = 1;
	int kt=0;
	int SomToDel,OtherSom;

	X=(int *)malloc((n+1)*sizeof(int));  
	L=(double *)malloc((n+1)*sizeof(double));
	Tree=(struct EDGE *)malloc((2*n-2)*sizeof(struct EDGE));
	Path=(struct EDGE *)malloc((n+2)*sizeof(struct EDGE));


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
	odp1(DI,X,&i,&j,n);

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

		if (DIS>2*epsilon)
		{
			while (S<DIS-(epsilon))
			{
				i=i+1;
				S=S+L[i];
			}
		}
		else { DIS=0; i=1; }

		Tree[p+1].U=n+k-1;
		Tree[p+1].V=Path[i].V;
		Tree[p+1].LN=S-DIS;
		if (Tree[p+1].LN<epsilon) Tree[p+1].LN=2*epsilon; 

		for (j=i+1;j<=P;j++)
		{
			Tree[p+j-i+1].U=Path[j].U;
			Tree[p+j-i+1].V=Path[j].V;
			Tree[p+j-i+1].LN=L[j];
			if (L[j]<2*epsilon) L[j]=2*epsilon; 
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
		if (fabs(Tree[i].LN-epsilon)<=2*epsilon)
			Tree[i].LN=0.0;
		ARETE[2*i-2]=Tree[i].U;
		ARETE[2*i-1]=Tree[i].V;
		LONGUEUR[i-1]=Tree[i].LN;   
		if (LONGUEUR[i-1]<2*epsilon) LONGUEUR[i-1] = 2*epsilon;
	} 
	/*/== rajoute le 22 avril 2005*/
	while((pasfini==1)&&(binaire==0)){
		//printf("non-binaire");
		//printf("\n");
		//for(i=1;i<=2*n-3-kt;i++){
		//	printf("\n%d : %d--%d -> %lf",i,ARETE[2*i-2],ARETE[2*i-1],LONGUEUR[i-1]);
		//}
		pasfini = 0;
		for (i=1;i<=2*n-3-kt;i++){
			if((LONGUEUR[i-1] == 2*epsilon)&&(ARETE[2*i-2]>n)&&(ARETE[2*i-1]>n)){	//= branche interne de taille=2*epsilon
				if(ARETE[2*i-2] > ARETE[2*i-1]){
					SomToDel = ARETE[2*i-2];
					OtherSom = ARETE[2*i-1];
				}
				else{
					SomToDel = ARETE[2*i-1];
					OtherSom = ARETE[2*i-2];
				}
				pasfini=1;

				/*/printf("\n\n edge to delete = %d \n\n",SomToDel);
				//== trouver les branches connexes*/
				for (j=1;j<=2*n-3-kt;j++){

					if(j!=i){
						if((ARETE[2*j-2] == SomToDel)){/*&&(ARETE[2*j-1] != ARETE[2*i-1])){*/
							ARETE[2*j-2] = OtherSom;
							/*/break;*/
						}
						else if((ARETE[2*j-1] == SomToDel)){/*/&&(ARETE[2*j-2] != ARETE[2*i-1])){*/
							ARETE[2*j-1] = OtherSom;
							/*/break;*/
						}
					}
				}
				for(j=i;j<=2*n-3-kt;j++){
					LONGUEUR[j-1] = LONGUEUR[j];
					ARETE[2*j-1] = ARETE[2*(j+1)-1];
					ARETE[2*j-2] = ARETE[2*(j+1)-2];
				}
				for(j=1;j<=2*n-3-kt-1;j++){
					if(ARETE[2*j-2] > SomToDel)
						ARETE[2*j-2]--;
					if(ARETE[2*j-1] > SomToDel)
						ARETE[2*j-1]--;
				}
				break;
			}
		} 
		if(pasfini)
			kt++;
	}

	/*printf("\n\n");
	for(i=1;i<=2*n-3-kt;i++){
	printf("\n%d : %d--%d -> %lf",i,ARETE[2*i-2],ARETE[2*i-1],LONGUEUR[i-1]);
	}
	*/


	/*printf("\n\n");*/

	free(X);
	free(Tree);
	free(L);
	free(Path);


	for (i=0;i<=n;i++)    
		free(D[i]);

	free(D);

	return kt;
}

//==============================================================
//
//==============================================================
int Bipartition_Table (double **D, int **B, int *PLACE, int n)
{
	int i,j,k,l,l1,*MaxCol,*X,EdgeNumberPath,m,uv,PlaceNumber,edge,*Path,M,F;
	double S,DIS,DIS1,*LengthPath;
	double EPS=1.e-5;
	double EPS1=1.e-2;

	/* Memory allocation */

	MaxCol=(int *)malloc((2*n)*sizeof(int));
	X=(int *)malloc((2*n+1)*sizeof(int));
	LengthPath=(double *)malloc((2*n)*sizeof(double));
	Path=(int *)malloc((2*n)*sizeof(int));

	/* Computation of a circular order X for D */     

	i=1; j=n; odp1(D,X,&i,&j,n);

	/* Initialization */ 
	for (i=1; i<=2*n-3; i++)
	{
		MaxCol[i]=0;
		PLACE[i]=0;
		for (j=1;j<=n;j++)
			B[i][j]=0;  
	}
	B[1][X[2]]=1; MaxCol[1]=X[2]; Path[1]=1; PlaceNumber=1;
	PLACE[1]=1; LengthPath[1]=D[X[1]][X[2]]; EdgeNumberPath=1; m=1;

	/* The main loop */

	for(k=2;k<=n-1;k++)
	{  
		/* Point 2.1 of the algorithm (see the referenced article by Makarenkov and Leclerc) */

		DIS=(D[X[1]][X[k]]+D[X[k]][X[k+1]]-D[X[1]][X[k+1]])/2; 
		DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
		
		//printf("\n\n%d,%d,%d\n\n",X[1],X[k],X[k+1]);
		if ((DIS<=-EPS1)||(DIS1<=-EPS1)) { printf("\n This is not an additive distance \n");
		free(MaxCol);free(X);free(LengthPath);free(Path);exit(1);return 0; }
		if (DIS<=EPS) DIS=0.0; if (DIS1<=EPS) DIS1=0.0;

		S=0.0; i=EdgeNumberPath; if (LengthPath[i]==0.0) i--;
		while (S<=DIS-EPS)
		{
			if (i==0) { S=DIS; break; }  /* checking the limit */  
			S=S+LengthPath[i];
			i--;
		}

		/* Point 2.2 of the algorithm */

		if (fabs(S-DIS)<=EPS) 
		{ 
			M=m+2; DIS=S;      
			if (i==0) F=1; 
			else if (i==EdgeNumberPath) F=2;
			else { M--; F=3; }   
		}    
		else {M=m+2; F=0;}


		if (M==m+2)
		{
			if (F==0) { uv=Path[i+1]; EdgeNumberPath=i+2; LengthPath[i+1]=S-DIS; LengthPath[i+2]=DIS1; 
			Path[i+1]=m+2; Path[i+2]=m+1;}
			else if (F==1) { uv=Path[1]; EdgeNumberPath=2; LengthPath[1]=0.0; LengthPath[2]=DIS1; 
			Path[1]=m+2; Path[2]=m+1;}
			else if (F==2) { uv=Path[EdgeNumberPath]; EdgeNumberPath=EdgeNumberPath+1;LengthPath[EdgeNumberPath]=DIS1; 
			Path[EdgeNumberPath-1]=m+2; Path[EdgeNumberPath]=m+1; } 

			for (j=1;j<=n;j++)
				B[m+2][j]=B[uv][j];   
			MaxCol[m+2]=MaxCol[uv];      
		}

		else 
		{
			EdgeNumberPath=i+1; LengthPath[i+1]=DIS1; Path[i+1]=m+1;
		}

		/* Point 2.3 of the algorithm */

		for (j=1;j<=EdgeNumberPath;j++)
			B[Path[j]][X[k+1]]=1;

		/* Point 2.4 of the algorithm */

		for (j=1;j<=EdgeNumberPath;j++)
			if (MaxCol[Path[j]]<X[k+1]) MaxCol[Path[j]]=X[k+1];  

		/* Point 2.5 of the algorithm */

		for (j=PlaceNumber;j>=1;j--) 
			PLACE[j+1]=PLACE[j];
		PLACE[1]=m+1; PlaceNumber++;

		if (M==m+2) {
			i=2; 
			while (PLACE[i]!=uv)
				i++;          
			for (j=PlaceNumber;j>=i+1;j--) 
				PLACE[j+1]=PLACE[j];
			PLACE[i+1]=m+2; PlaceNumber++;}

		i=M-1; edge=2;
		do 
		{
			if (PLACE[i]==Path[edge]) 
			{
				edge++; j=i+1; 
				while (X[k+1]>MaxCol[PLACE[j]]) 
					j++; 
				if (j>i+1) 
				{            
					l1=PLACE[i];
					for (l=i+1;l<=j-1;l++) 
						PLACE[l-1]=PLACE[l];              
					PLACE[j-1]=l1;

				} 
			}
			i--;
		} while (i!=0); 

		m=M;
	}


	/* memeory liberation */

	free(MaxCol);
	free(X);
	free(LengthPath);
	free(Path);

	return m;

}

//==============================================================
//
//==============================================================
int Table_Comparaison (int **B, int ** B1, int *PLACE, int *PLACE1, int m, int m1,int n)
{
	int RF=0,i,p,p1;

	p=1; p1=1;

	while ((p<=m)&&(p1<=m1))
	{
		i=n;
		while ((B[PLACE[p]][i]==B1[PLACE1[p1]][i])&&(i>1))
			i--;
		if (i==1) { RF=RF+1; p++; p1++; }
		else if (B[PLACE[p]][i]>B1[PLACE1[p1]][i]) p1++;
		else p++; 

	}
	RF=(m-RF)+(m1-RF);

	return RF;
}


void approx_arb(double **DISS,double **TM,double **TMnew,double **W, int *Iternum, long int *ARETE, double *LONGUEUR,int choix,int *kt,int binaire,int n)
{

	long int  *Level, *Score, *EX1, *EX2, *Succ1, *Succ2, *Succ11, *Succ22, *Neighbour1, *Neighbour2;
	int i,j,k,k1,j1,i1,i2,*Flag, **Part, **Vertices, *VertexNumber;
	long int *ARETE1;                                       
	int *degree=0, a=0, exit_number;    
	double *L, **B, *C, Sum, Sum1, *Path, EQ, l;
	int na1,ns1,na,ns;

	/*/ Variable declaration*/
	//printf("\n========== %d",(*kt));;
	Level=(long int *)malloc((2*n-1)*sizeof(long int));  
	Score=(long int *)malloc((2*n-1)*sizeof(long int));  
	EX1=(long int *)malloc((2*n-1)*sizeof(long int));  
	EX2=(long int *)malloc((2*n-1)*sizeof(long int));
	Succ1=(long int *)malloc((2*n-1)*sizeof(long int)); 
	Succ2=(long int *)malloc((2*n-1)*sizeof(long int));
	Succ11=(long int *)malloc((2*n-1)*sizeof(long int)); 
	Succ22=(long int *)malloc((2*n-1)*sizeof(long int));
	/*/ARETE=(long int *)malloc((4*n-2)*sizeof(long int));*/
	C=(double *)malloc((2*n-1)*sizeof(double));
	B=(double **)malloc((2*n-1)*sizeof(double*));

	VertexNumber=(int *)malloc((2*n-1)*sizeof(int));    
	Vertices=(int **)malloc((2*n-1)*sizeof(int*));
	Part=(int **)malloc((2*n-1)*sizeof(int*));

	Neighbour1=(long int *)malloc((2*n-1)*sizeof(long int));
	Neighbour2=(long int *)malloc((2*n-1)*sizeof(long int));
	Flag=(int *)malloc((2*n-1)*sizeof(int));

	L=(double *)malloc((2*n-1)*sizeof(double));
	/*/LONGUEUR=(double *)malloc((2*n-1)*sizeof(double));*/
	Path=(double *)malloc((n+1)*sizeof(double));  

	for (i=0;i<2*n-1;i++){
		Score[i]=0;	
	}
	for (i=0;i<=2*n-2;i++)
	{
		Vertices[i]=(int*)malloc((n+1)*sizeof(int));  
		B[i]=(double*)malloc((2*n-1)*sizeof(double));
		Part[i]=(int*)malloc((2*n-1)*sizeof(int));   
		if ((B[i]==NULL)||(Part[i]==NULL)||(Vertices[i]==NULL))
		{
			printf("Data matrix is too large"); 
			return;
		}
	}   

	/*/ Variable initialisation*/

	if(choix==0){
		*kt = Tree_edges (TM,ARETE,LONGUEUR,n,binaire);
		for(i=0;i<=n;i++)
			for(j=0;j<=n;j++)
				TMnew[i][j] = TM[i][j]; 
	}
	else{
		for(i=0;i<=2*n-2;i++)
			for(j=0;j<=2*n-2;j++)
				TMnew[i][j] = TM[i][j]; 
	}

	/*/return; */
	na1 = 2*n-3-*kt;
	ns1 = 2*n-2-*kt; 


	na = 2*n-3;   /*/start of the new block - june 2005*/
	ns = 2*n-2;                                         
	for (i=2*n-3-*kt+1; i<=2*n-3; i++)                      
		LONGUEUR[i-1] = 5*epsilon;  

	ARETE1	=(long int*)malloc(2*(2*n-2)*sizeof(long int)); 
	degree	=(int*)malloc((2*n-1)*sizeof(int)); 

	for (i=1; i<=na1; i++)                                  
	{        
		ARETE1[2*i-2] = ARETE[2*i-2];                         
		ARETE1[2*i-1] = ARETE[2*i-1];                         
	}

	i=1;
	while (i<=*kt)
	{
		for (j=1; j<=ns; j++)  //-                               
			degree[j] = 0;

		for (j=1; j<=na1+i-1; j++)    
		{
			if (degree[ARETE1[2*j-2]]==3) { exit_number = ARETE1[2*j-2]; break; }
			else ++degree[ARETE1[2*j-2]];

			if (degree[ARETE1[2*j-1]]==3) { exit_number = ARETE1[2*j-1]; break; }
			else ++degree[ARETE1[2*j-1]];
		}

		a=0;
		for (j=1; j<=na1+i-1; j++) 
		{
			if (ARETE1[2*j-2]==exit_number) { ARETE1[2*j-2]=ns1+i; a++; }
			if (ARETE1[2*j-1]==exit_number) { ARETE1[2*j-1]=ns1+i; a++; }
			if (a==2) break;
		}

		ARETE1[2*(na1+i)-2] = exit_number;
		ARETE1[2*(na1+i)-1] = ns1+i;
		i++;
	}   /*/end of the new block - june 2005      */                                                        


	for (i=1;i<=na;i++) //-
	{        
		EX1[i]=ARETE1[2*i-2]; /*/modified line - june 2005*/
		EX2[i]=ARETE1[2*i-1]; /*/modified line - june 2005*/

		L[i]=LONGUEUR[i-1];

		Flag[i]=0;

		if (EX1[i]<=n)
		{ Score[EX1[i]]=3; Succ1[EX1[i]]=0; Succ2[EX1[i]]=0; Score[EX2[i]]=0; }

		else if (EX2[i]<=n)
		{ Score[EX2[i]]=3; Succ1[EX2[i]]=0; Succ2[EX2[i]]=0; Score[EX1[i]]=0; }

		else { Score[EX1[i]]=0; Score[EX2[i]]=0; } 
	}

	for (i=1;i<=na;i++)//-
	{ 
		if (i<=n) { Path[i]=0.0; TMnew[i][i]=0; }
		if ((EX1[i]<=n)||(EX2[i]<=n))
		{
			VertexNumber[i]=1;
			if (EX1[i]<=n) Vertices[i][1]=EX1[i];
			else Vertices[i][1]=EX2[i];     
		}
		else VertexNumber[i]=0;

		for (j=1;j<=na;j++)//-
		{
			Part[i][j]=0;
			Part[j][i]=0;   
		} 
	}


	/*/   Filling in of the vectors Level(2n-3), VertexNumber(2n-3),
	//   and the matrices Part(2n-3, 2n-3), Vertices(2n-3,n)   */ 

	int cpt=0;
	j=1; 
	while (j<ns)
	{    
		cpt++;
		//printf("%d ",cpt);
		if(cpt > 10000) {printf("\nProbleme approx_arb"); exit(-1);}
		for (i=1; i<=na; i++)//-
		{
			if (Score[EX1[i]]==5) 
				Score[EX1[i]]=3;
			if (Score[EX2[i]]==5) 
				Score[EX2[i]]=3;

			if (((Score[EX1[i]]==3)&&(Score[EX2[i]]!=4))||((Score[EX2[i]]==3)&&(Score[EX1[i]]!=4))) 
			{ 
				if (Score[EX1[i]]==3) 
				{ k=EX2[i]; k1=EX1[i]; }
				else 
				{ k=EX1[i]; k1=EX2[i]; }

				if (Score[k]<2) 
				{
					Score[k]=Score[k]+1;
					if (Score[k]==1) 
					{ Succ1[k]=i; Neighbour1[k]=k1; Neighbour2[k]=0;}
					else  
					{ Succ2[k]=i; Neighbour2[k]=k1; }
				}
			}
		}

		for (i=1; i<=na; i++)//-
		{
			/*/if(j>ns) break;*/
			if (((Score[EX1[i]]==1)&&(Score[EX2[i]]==3))||((Score[EX1[i]]==2)&&(Score[EX2[i]]==3))||
				((Score[EX1[i]]==3)&&(Score[EX2[i]]==1))||((Score[EX1[i]]==3)&&(Score[EX2[i]]==2)))

			{ 

				if ((Score[EX1[i]]==1)||(Score[EX1[i]]==2)) { k=EX1[i]; k1=EX2[i]; }
				else { k=EX2[i]; k1=EX1[i]; } 

				if(Flag[i]==0) 
				{ 
					Succ11[i]=Succ1[k1];
					Succ22[i]=Succ2[k1];
					Level[j]=i; j=j+1; Score[k1]=4;
					Flag[i]=1;
		//			printf("\nk=%d , Neighbour1[k] = %d , Neighbour2[k] = %d , %d-%d", k,Neighbour1[k],Neighbour2[k],Score[Neighbour1[k]],Score[Neighbour2[k]]);
					if ((Score[Neighbour1[k]]==4)&&(Score[Neighbour2[k]]==4)) Score[k]=5;

					if (j>n+1)
					{              
						VertexNumber[i]=VertexNumber[Succ11[i]]+VertexNumber[Succ22[i]];
						for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)
							Vertices[i][i1]=Vertices[Succ11[i]][i1];

						for (i1=VertexNumber[Succ11[i]]+1; i1<=VertexNumber[i]; i1++)
							Vertices[i][i1]=Vertices[Succ22[i]][i1-VertexNumber[Succ11[i]]];


						Part[Succ11[i]][i]=1; Part[Succ22[i]][i]=1;
						Part[i][Succ11[i]]=1; Part[i][Succ22[i]]=1; 
						for (i1=1; i1<=na; i1++)
						{          
							if ((Part[Succ11[i]][i1]==1)||(Part[Succ22[i]][i1]==1))
							{
								Part[i][i1]=1; Part[i1][i]=1;
							}
						}
					}          
				}

				/*/ if(j>ns) break;*/

				if(Score[k]==2)
				{ 
					Score[k]=5;
					Succ11[Succ2[k]]=Succ1[Neighbour2[k]];
					Succ22[Succ2[k]]=Succ2[Neighbour2[k]];
					Level[j]=Succ2[k]; j=j+1; Score[Neighbour2[k]]=4;
					Flag[Succ2[k]]=1;       

					j1=i;
					i=Succ2[k];

					if (j>n+1)
					{              
						VertexNumber[i]=VertexNumber[Succ11[i]]+VertexNumber[Succ22[i]];
						for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)
							Vertices[i][i1]=Vertices[Succ11[i]][i1];

						for (i1=VertexNumber[Succ11[i]]+1; i1<=VertexNumber[i]; i1++)
							Vertices[i][i1]=Vertices[Succ22[i]][i1-VertexNumber[Succ11[i]]];


						Part[Succ11[i]][i]=1; Part[Succ22[i]][i]=1;
						Part[i][Succ11[i]]=1; Part[i][Succ22[i]]=1; 
						for (i1=1; i1<=na; i1++)
						{          
							if ((Part[Succ11[i]][i1]==1)||(Part[Succ22[i]][i1]==1))
							{
								Part[i][i1]=1; Part[i1][i]=1;
							}
						}
					} 
					i=j1;            

				} 
			} 
		} 

		/*/if(j>ns) break;*/

		for (i=1; i<=na; i++)//-
		{ 
			/*/  if(j>ns) break;*/

			if (((Score[EX1[i]]==5)&&(Score[EX2[i]]==3))||((Score[EX1[i]]==3)&&(Score[EX2[i]]==5))||
				((Score[EX1[i]]==5)&&(Score[EX2[i]]==5)))
			{      
				if ((Score[EX1[i]]==3)||((Score[EX1[i]]==5)&&(Score[EX2[i]]==5)))
				{ k=EX1[i]; k1=EX2[i]; }
				else { k=EX2[i]; k1=EX1[i]; } 
				Level[j]=i; j=j+1;    
				Succ11[i]=Succ1[k];
				Succ22[i]=Succ2[k];

				Score[k]=4;
				Score[k1]=4;

				if (j>n+1)
				{   
					//printf("\n============ %d ",i);           
					//printf("\n============ %d , %d , %d ",i,Succ11[i],Succ22[i]);           
					if(Succ22[i] > INFINI || Succ11[i] > INFINI){
						printf("\nErreur dans approx_arb (utils_tree.cpp:1111)\n");
						exit(-1);
					}
					VertexNumber[i]=VertexNumber[Succ11[i]]+VertexNumber[Succ22[i]];
					for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)
						Vertices[i][i1]=Vertices[Succ11[i]][i1];

					for (i1=VertexNumber[Succ11[i]]+1; i1<=VertexNumber[i]; i1++)
						Vertices[i][i1]=Vertices[Succ22[i]][i1-VertexNumber[Succ11[i]]];
					Part[Succ11[i]][i]=1; Part[Succ22[i]][i]=1;
					Part[i][Succ11[i]]=1; Part[i][Succ22[i]]=1; 
					for (i1=1; i1<=na; i1++)
					{          
						if ((Part[Succ11[i]][i1]==1)||(Part[Succ22[i]][i1]==1))
						{
							Part[i][i1]=1; Part[i1][i]=1;
						}
					}
				}            
			}       
		}
	}  


	/*/ Filling in of the matrix B(2n-3,2n-3) and the vector C(2n-3)  */

	for (i=1; i<=na; i++)//-
	{
		B[i][i]=0;
		C[i]=0;
		for (j=1; j<=na; j++)//-
		{   
			if ((EX1[i]<=n)||(EX2[i]<=n))
			{
				if (EX1[i]<=n) k=EX1[i];
				else k=EX2[i];
				if (((EX1[j]<=n)||(EX2[j]<=n))&&(i!=j))       
				{           
					if (EX1[j]<=n) k1=EX1[j];
					else k1=EX2[j];

					B[i][j]=1.0; //W[k][k1];
					B[j][i]=1.0; //W[k][k1];
					B[i][i]=B[i][i]+B[i][j];
					C[i]=C[i]+DISS[k][k1]; //*W[k][k1];     
				}
			}
		}
	}

	for (k=n+1; k<=na; k++)//-
	{   
		i=Level[k];
		Sum=0.0;
		Sum1=0.0;
		for (j=1; j<=VertexNumber[Succ11[i]]; j++) 
		{
			for (j1=1; j1<=VertexNumber[Succ22[i]]; j1++)     
			{
				Sum=Sum+2*1.0;//W[Vertices[Succ11[i]][j]][Vertices[Succ22[i]][j1]];
				Sum1=Sum1+2*DISS[Vertices[Succ11[i]][j]][Vertices[Succ22[i]][j1]]*1.0;
					//W[Vertices[Succ11[i]][j]][Vertices[Succ22[i]][j1]];
			}
		}     
		B[i][i]=B[Succ11[i]][Succ11[i]]+B[Succ22[i]][Succ22[i]]-Sum;
		C[i]=C[Succ11[i]]+C[Succ22[i]]-Sum1;  
	}


	for (j1=n+1; j1<=na; j1++)//-
	{   
		i=Level[j1];
		for (i1=1; i1<j1; i1++)
		{       
			j=Level[i1];
			if ((j==Succ11[i])||(j==Succ22[i]))
			{  
				if (j==Succ11[i]) k=Succ22[i];
				else k=Succ11[i]; 
				B[i][j]=(B[i][i]+B[j][j]-B[k][k])/2;
			}  

			else if (Part[i][j]==1)
			{        
				if (Part[j][Succ11[i]]==1)
					B[i][j]=B[j][Succ11[i]]-B[j][Succ22[i]];
				else 
					B[i][j]=B[j][Succ22[i]]-B[j][Succ11[i]];    
			}

			else    
			{        
				B[i][j]=B[Succ11[i]][j]+B[Succ22[i]][j];        
			}
			B[j][i]=B[i][j];

		}
	}  


	/* Gausse-Seidel iteretive procedure to polish vector of lengths L(2n-3)*/

	for(k=1;k<=*Iternum;k++)
	{    


		for(i=1;i<=na1;i++) /*/modified line - june 2005*/
		{
			EQ=0;
			Sum=0;
			Sum1=0;
			for (j=1;j<=na1;j++) /*/modified line - june 2005*/
			{
				if (j>=i+1) Sum=Sum+B[i][j]*L[j];
				if (j<=i-1) Sum1=Sum1+B[i][j]*L[j];
			}
			l=L[i];
			L[i]=(-Sum-Sum1+C[i])/B[i][i];
			if (L[i]<5*epsilon) L[i]=5*epsilon;
			EQ=EQ+sqrt((double)((L[i]-l)*(L[i]-l)));   
		}    
	}


	for (i=1;i<=na1;i++) /*/modified line - june 2005*/
	{
		LONGUEUR[i-1]=L[i];
		if (L[i]<5*epsilon) L[i] = LONGUEUR[i-1] = 5*epsilon;
	}


	/*/ Computing of new tree distance matrix TM(n,n) from the list of edges*/

	for (j=n+1; j<=na; j++)//-
	{   
		i=Level[j];

		for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)   
			Path[Vertices[Succ11[i]][i1]]=Path[Vertices[Succ11[i]][i1]]+L[Succ11[i]];  

		for (i1=1; i1<=VertexNumber[Succ22[i]]; i1++)     
			Path[Vertices[Succ22[i]][i1]]=Path[Vertices[Succ22[i]][i1]]+L[Succ22[i]];   

		for (i1=1; i1<=VertexNumber[Succ11[i]]; i1++)   
		{   
			k=Vertices[Succ11[i]][i1];
			for (j1=1; j1<=VertexNumber[Succ22[i]]; j1++)     
			{            
				k1=Vertices[Succ22[i]][j1];
				TMnew[k][k1]=Path[k]+Path[k1];
				TMnew[k1][k]=TMnew[k][k1];
			}
		}              
	}

	i=Level[na];
	if ((EX1[i]==EX1[Succ11[i]])||(EX1[i]==EX2[Succ11[i]])) k=EX2[i];
	else k=EX1[i];

	j1=1;
	for (i1=1; i1<=na-1; i1++)//-
	{    
		k1=Level[i1];
		if ((EX1[k1]==k)||(EX2[k1]==k))
		{
			if (j1==1) { i=k1; j1=2; } 
			else { j=k1; break; }  
		}
	}

	for (i1=1; i1<=VertexNumber[i]; i1++)   
	{    
		for (j1=1; j1<=VertexNumber[j]; j1++)     
		{   
			k=Vertices[i][i1];
			k1=Vertices[j][j1];
			//printf("\n%d-%d path=%lf",k,k1,Path[k]);
			TMnew[k][k1]=Path[k];
			
			TMnew[k][k1]+=Path[k1];
			TMnew[k][k1]+=L[i];
			TMnew[k][k1]+=L[j];
			TMnew[k1][k]=TMnew[k][k1];
		}
	}        

	i2=Level[na];//-
	for (i1=1; i1<=VertexNumber[i2]; i1++)   
	{             
		k=Vertices[i2][i1];
		for (j1=1; j1<=VertexNumber[i]; j1++)     
		{   
			k1=Vertices[i][j1];
			TMnew[k][k1]=Path[k]+Path[k1]+L[i2]+L[i];
			TMnew[k1][k]=TMnew[k][k1];
		}     
		for (j1=1; j1<=VertexNumber[j]; j1++)     
		{   
			k1=Vertices[j][j1];
			TMnew[k][k1]=Path[k]+Path[k1]+L[i2]+L[j];
			TMnew[k1][k]=TMnew[k][k1];
		}          
	}     

	free(Level);  
	free(Score);  
	free(EX1);  
	free(EX2);
	free(Succ1); 
	free(Succ2);
	free(Succ11); 
	free(Succ22);
	free(C);  
	free(VertexNumber);     
	free(Neighbour1);
	free(Neighbour2);
	free(Flag);  
	free(L);
	free(Path);

	for (i=0;i<=ns;i++)
	{ 
		free(B[i]);
		free(Part[i]);
		free(Vertices[i]);  	   
	} 
	free(B);
	free(Part);
	free(Vertices); 

	free(ARETE1); /*/new line june 2005*/
	free(degree); /*/new line june 2005*/
}

/*=======================================
//
//=======================================*/
int findFils(double ** Adjacence,int sommet,int n){

	int i;
	
	for(i=1;i<=2*n-2;i++){
		if(Adjacence[sommet][i]<INFINI){
			Adjacence[sommet][i] = Adjacence[i][sommet] = INFINI;
			return i;
		}
	}
	return -1;
}

/*=======================================
//
//=======================================*/
struct TNoeud * CreerSousArbre(long int *ARETE,int *indice, double ** Adjacence,int sommet,int n){

	int i=0;
	int *tableau;
	int nbElt=0;
	int nouveauFils;
	struct TNoeud * node;

	if(sommet==-1)
		return NULL;

	tableau = (int*)malloc(100*sizeof(int));

	node = (struct TNoeud*)malloc(sizeof(struct TNoeud));
	node->NoNoeud = sommet;
	node->nbfils = 0;
	node->fils = (struct TNoeud**)malloc(100*sizeof(struct TNoeud*));
	for(i=0;i<100;i++)
		node->fils[i] = NULL;

	nouveauFils = -1;

	do{
		nouveauFils = findFils(Adjacence,sommet,n);
		if(nouveauFils != -1){
			tableau[nbElt] = nouveauFils;
			nbElt++;
		}
	}while(nouveauFils != -1);

	for(i=0;i<nbElt;i++){	
		node->fils[node->nbfils] = CreerSousArbre(ARETE,indice,Adjacence,tableau[i],n);
		node->nbfils = node->nbfils + 1;
	}
	free(tableau);
	return node;
}

/*/=============================================
//= tri bulle d'un tableau d'entiers
//=============================================*/
void sortIntTab(int * tab, int debut, int fin){
	
	int i,j,tmp;

	for(i=debut;i<=fin;i++){
		for(j=i+1;j<=fin;j++){
			if(tab[i] > tab[j]){
				tmp = tab[i];
				tab[i] = tab[j];
				tab[j] = tmp;
			}
		}
	}
}

/*/============================================================
//
//============================================================*/
void ParcoursArbre(struct TNoeud * unNoeud, struct DescTree * SousMatriceTree){

	int j,k;
	int somme=0;
	int nbSommet=0;
	int nbfils;
	
	if(unNoeud != NULL){
		nbfils = unNoeud->nbfils;
		/*printf("\n%d",unNoeud->NoNoeud);*/


		if(nbfils != 0){

			for(j=0; j<nbfils;j++){
				ParcoursArbre(unNoeud->fils[j],SousMatriceTree);
			}
		
			for(j=0; j<nbfils;j++){
				/*toto = unNoeud->fils[j]->NoNoeud;*/
				somme = somme + SousMatriceTree[unNoeud->fils[j]->NoNoeud].nbSommet;
			}
			SousMatriceTree[unNoeud->NoNoeud].nbSommet = somme; /*SousMatriceTree[unNoeud->droit->NoNoeud].nbSommet + SousMatriceTree[unNoeud->gauche->NoNoeud].nbSommet;*/
			SousMatriceTree[unNoeud->NoNoeud].Tableau = (int*)malloc((somme+1)*sizeof(int));
			
			for(j=0;j<nbfils;j++){
				for(k=1;k<=SousMatriceTree[unNoeud->fils[j]->NoNoeud].nbSommet;k++) {
					nbSommet++;	
					SousMatriceTree[unNoeud->NoNoeud].Tableau[nbSommet] = SousMatriceTree[unNoeud->fils[j]->NoNoeud].Tableau[k];
				}
			}
			/*sort(listeSommet.begin(),listeSommet.end());*/
			sortIntTab(SousMatriceTree[unNoeud->NoNoeud].Tableau,1,nbSommet);

			/*SousMatriceTree[unNoeud->NoNoeud].nbSommet = listeSommet.size();*/
			SousMatriceTree[unNoeud->NoNoeud].nbSommet = nbSommet;
			
		}
		else{

			SousMatriceTree[unNoeud->NoNoeud].Tableau = (int*)malloc((2)*sizeof(int));

			SousMatriceTree[unNoeud->NoNoeud].Tableau[1] = unNoeud->NoNoeud;
			SousMatriceTree[unNoeud->NoNoeud].nbSommet = 1;
		}
	}
	
}

/*=======================================
//
//=======================================*/
void viderArbre(struct TNoeud * A){
	int i;
	
	if(A->nbfils != 0){
		for(i=0; i<A->nbfils;i++){
			viderArbre(A->fils[i]);
		}
	//	printf("%d ",A->NoNoeud);
		free(A->fils);
		free(A);
	}
	else{
		free(A->fils);
		free(A);
	}
}

/*=======================================
//
//=======================================*/
void AfficherArbre(struct TNoeud * A,int prof){
	
	int i,j, taille;
	
	if(A != NULL){
		if(A->nbfils== 0){
			for(i=0;i<prof;i++) printf("  ");
			printf("-- %d\n",A->NoNoeud);
		}
		else{
			taille = A->nbfils;
			for(i=0; i<taille/2;i++){
				AfficherArbre(A->fils[i],prof+2);
				printf("\n");
				//for(j=0; j<prof;j++) printf(" ");
				//printf("|");
			}
			for(i=0;i<prof;i++) printf("  ");
			printf("-- %d\n",A->NoNoeud);
			for(i=taille/2; i<taille;i++){
				AfficherArbre(A->fils[i],prof+2);
			}
		}
	}
}

/*=======================================
//
//=======================================*/
void AfficherArbre2(struct TNoeud * A,int prof,char ** NOMS){
	
	int i,j, taille;
	
	if(A != NULL){
		if(A->nbfils== 0){
			for(i=0;i<prof;i++) printf("  ");
			printf("-- %s\n",NOMS[A->NoNoeud]);
		}
		else{
			taille = A->nbfils;
			for(i=0; i<taille/2;i++){
				AfficherArbre2(A->fils[i],prof+2,NOMS);
				printf("\n");
			}
			for(i=0;i<prof;i++) printf("  ");
			printf("-- %s\n",NOMS[A->NoNoeud]);
			for(i=taille/2; i<taille;i++){
				AfficherArbre2(A->fils[i],prof+2,NOMS);
			}
		}
	}
}

/*/=================================================================================================
//
//=================================================================================================*/
void printTree(long int *ARETE,double ** DIST,int P,double ** Adjacence2,struct DescTree *DT,int n,int kt,char ** NOMS)
{
	/*/===== variables =====*/
	int i,j,k,l,taille,plus,ns=2*n-2;
	double ** Adjacence =(double**)malloc((2*n+1)*sizeof(double*));
	struct TNoeud * arbre;
	
	for (i=0;i<=2*n;i++)
		Adjacence[i]=(double*)malloc((2*n+1)*sizeof(double));

	/*/==== Processus =====*/

	for(i=0;i<=2*n-2;i++)
		for(j=0;j<=2*n-2;j++)
			Adjacence[i][j] = Adjacence2[i][j];
	i=1;

	Adjacence[P][n] = Adjacence[n][P] = INFINI; 

	arbre = CreerSousArbre(ARETE,&i,Adjacence,P,n);

	printf("\n-------------------------------\n");
	AfficherArbre2(arbre,0,NOMS);
	printf("\n-------------------------------");
}

/*/=================================================================================================
//
//=================================================================================================*/
void RechercherBipartition(long int *ARETE,double ** DIST,int P,double ** Adjacence2,struct DescTree *DT,int n,int kt)
{
	/*/===== variables =====*/
	int i,j,k,l,taille,plus,ns=2*n-2;
	double ** Adjacence =(double**)malloc((2*n+1)*sizeof(double*));
	struct TNoeud * arbre;
	
	for (i=0;i<=2*n;i++)
		Adjacence[i]=(double*)malloc((2*n+1)*sizeof(double));

	/*/==== Processus =====*/

	for(i=0;i<=2*n-2;i++)
		for(j=0;j<=2*n-2;j++)
			Adjacence[i][j] = Adjacence2[i][j];
	i=1;

	Adjacence[P][n] = Adjacence[n][P] = INFINI; 

	arbre = CreerSousArbre(ARETE,&i,Adjacence,P,n);

	printf("\n-------------------------------\n");
	AfficherArbre(arbre,0);
	printf("\n-------------------------------");


	ParcoursArbre(arbre,DT);
	//printf("\nRechercher bipartition :");
	for(i=1;i<=ns-kt;i++){
    	if(i!=n){
        	taille = DT[i].nbSommet;

		/*	printf("\n=== %d\n",i);
			for(j=1;j<=DT[i].nbSommet;j++)
				printf(" %d",DT[i].Tableau[j]);*/

        	/*if(taille == 3) plus=1;
                else */
					plus=1;
			
            DT[i].Matrice = (double**)malloc((taille+1+plus)*sizeof(double*));
            for(j=0;j<=taille+plus;j++){
                DT[i].Matrice[j] = (double*)malloc((taille+1+plus)*sizeof(double));
            }
            for(k=1;k<=taille;k++){
                for(l=1;l<=taille;l++){
                    DT[i].Matrice[k][l] = DIST[DT[i].Tableau[k]][DT[i].Tableau[l]] ;
                }
            }
			 
			for(l=1;l<=taille;l++){
				DT[i].Matrice[taille+1][l] = DT[i].Matrice[l][taille+1] = DIST[n][DT[i].Tableau[l]];
            }

            DT[i].Matrice[taille+1][taille+1] = 0.0;
          /*  if(taille == 3){
                for(l=1;l<=taille;l++){
                    DT[i].Matrice[taille+1][l] = DT[i].Matrice[l][taille+1] = DIST[n][DT[i].Tableau[l]];
                }
                DT[i].Matrice[taille+1][taille+1] = 0.0;
            }*/
     	}
	}
	//printf("\n Fin Rechercher bipartition :");
	viderArbre(arbre);
	for (i=0;i<=2*n;i++)
		free(Adjacence[i]);
	free(Adjacence);
}

static void xtoa (unsigned long val,char *buf,unsigned radix,int is_neg){
	char *p;                /* pointer to traverse string */
	char *firstdig;         /* pointer to first digit */
	char temp;              /* temp char */
	unsigned digval;        /* value of digit */

	p = buf;

	if (is_neg) {
		/* negative, so output '-' and negate */
		*p++ = '-';
		val = (unsigned long)(-(long)val);
	}

	firstdig = p;           /* save pointer to first digit */

	do {
		digval = (unsigned) (val % radix);
		val /= radix;       /* get next digit */

		/* convert to ascii and store */
		if (digval > 9)
			*p++ = (char) (digval - 10 + 'a');  /* a letter */
		else
			*p++ = (char) (digval + '0');       /* a digit */
	} while (val > 0);

	/* We now have the digit of the number in the buffer, but in reverse
	order.  Thus we reverse them now. */

	*p-- = '\0';            /* terminate string; p points to last digit */

	do {
		temp = *p;
		*p = *firstdig;
		*firstdig = temp;   /* swap *p and *firstdig */
		--p;
		++firstdig;         /* advance to next two digits */
	} while (firstdig < p); /* repeat until halfway */
}

/* Actual functions just call conversion helper with neg flag set correctly,
and return pointer to buffer. */

char * itoa_(int val,char *buf,int radix){
	if (radix == 10 && val < 0)
		xtoa((unsigned long)val, buf, radix, 1);
	else
		xtoa((unsigned long)(unsigned int)val, buf, radix, 0);
	return buf;
}


//===============================================================================================
//
//===============================================================================================
int lectureNewick(string newick, long int * ARETE, double * LONGUEUR, char ** lesNoms, int *kt)
{
	int n;                                     
	int cpt_x;
	
	// Ce sous programme permet de lire un arbre au format newick et de le transcrire dans
	// des vecteurs arete-longueur commencant � 1
	// ATTENTION: les noms commencent � 0
	// 7 octobre 2004
	// Elmaestro

	// TODO: Add your command handler code here
	int i,j,j1,k, a, a1, a2,a3, VertexNumber,numero;
	char symbol, *string, *string1, *string2, *string4;
	int taxaPos; // le nombre de taxas recup�r�
	int aretePos; // le nombre d'aretes recup�r�
	char symbolOld =' ';
	int zz, xx,jj,ii;
	double longueur;
	char * tempString;
	int cpt=0;
	//char *string4 = (char*) malloc(100000 * sizeof(char));
	string4 = (char*)malloc((100000) * sizeof(char));
	//char *string4 = (char*) malloc(100000);
	//char *string4 = new char[10000];

	int temoin=0;
	int cpt_parenthese=0;
	//long int *ARETE;
	//double *LONGUEUR; 	

	//Correctness of the Newick format verification
	i=0;
	n = 0;

	do
	{
		symbol = newick.at(cpt);
		cpt++;
		if (symbol==':') i++;
		if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}  while(cpt < newick.size());

	cpt=0;
	if(i<=2*n-3)(*kt)=i;
	else (*kt)=2*n-3;

	if (i == 0) { printf("Incorrect Newick file format. Edge lengths must be indicated after a ':' characters."); exit(FAIL);}

	if ((i>2*n-3) || (i<n)){
		//printf("Unable to read your input data, please check it and try again...");
		//exit(-1);
	}

	j=0;
	do{ 
	   symbol=newick.at(cpt); 
	   cpt++;
	   if (symbol=='(') j++; 
	}while(cpt < newick.size());
	
	cpt=0;	
	j1=0;
	
	do{
		symbol=newick.at(cpt);
		cpt++;
		if (symbol==')') j1++; 
	}while(cpt < newick.size());  
	
	cpt=0;

	// verification des arit�s de l'arbre
	if (j1!=j) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); exit(FAIL);}
	//else if (j!=n-2) { printf("Incorrect Newick file format. Only trees with vertices of degree 1 and 3 are allowed by T-REX."); fclose (data); exit(FAIL);}

	k=0;

	do{ 
		symbol=newick.at(cpt);
		cpt++;
		if (symbol==',') k++; 
	}while(cpt < newick.size());   
	
	cpt=0;
	//if (k!=(n-1)) { printf("Incorrect Newick file format. Number of objects must be equal to number of commas plus 1."); fclose (data); exit(FAIL);}

	a=0;

	do{ 
		symbol=newick.at(cpt);
		cpt++;
		if (symbol==';') a++; 
	}while(cpt < newick.size());    
	cpt=0;

	if (a==0) { printf("Incorrect Newick file format. Newick string must be followed by a ';' character."); exit(FAIL);}
	else if (a>1) { printf("Incorrect Newick file format. Newick string must contain (in the end) only one ';' character."); exit(FAIL);}

	a=0;
	do{
		symbol=newick.at(cpt);
		cpt++;
		a++;
	}while(symbol == ' '); 
	
	cpt=0;

	if (symbol!='(') { printf("Incorrect Newick file format. Newick string must begin with a '(' character."); exit(FAIL);}  

	a=0;
	
	do{ 
		symbol=newick.at(cpt);
		cpt++;
		if (symbol=='%') a++; 
	}while(cpt < newick.size());

	cpt=0;

	if (a>0) { printf("Incorrect Newick file format. Newick string cannot contain \'%%\' character."); exit(FAIL);}

	do
	{
		symbol=newick.at(cpt);
		cpt++;
		if ((symbol=='(')||(symbol==','))		
		{ 
			symbol=newick.at(cpt);
			cpt++;
			a=0;
			if ((symbol!='(')&&(symbol!=',')&&(symbol!=';')&&(symbol!=':'))		   
			{ 
				cpt--; 
				do{
					symbol=newick.at(cpt);
					cpt++;
					a++;
				}while ((symbol!=':')&&(cpt < newick.size()));  
			}
			else cpt--;
			if (a>50) { printf("Incorrect Newick file format. Names of objects must not exceed 50 characters.");  exit(FAIL);}   
		}	
	}while(cpt < newick.size());
	cpt=0;

	string = (char*)malloc((100000) * sizeof(char));
	string2 = (char*)malloc((100000) * sizeof(char));
	//string3 = (char*)malloc((1000*n) * sizeof(char));
	string1 = (char*)malloc((100000) * sizeof(char));

	for(int cpt_string=0; cpt_string<100000; cpt_string++){
		string[cpt_string] = ' ';
		string2[cpt_string] = ' ';
	}

	if ((string == NULL)||(string1 == NULL)||(string2 == NULL)/*||(string3 == NULL)*/)
	{ printf("Input data are too large or not a correct Newick file chosen"); exit(FAIL);}

	a=0;

	do{		
		symbol=newick.at(cpt);
		cpt++;
		if ((symbol!=' ')&&(symbol!='\n')&&(symbol!='\t')) { string[a++]=symbol; } 
	}while(cpt < newick.size()); 

	int boot_value;
	int temoin3 = 0;
	k=0; VertexNumber=n;
	//a1 = 0;
	//a2 = 0;
	taxaPos =1;    // nous allons commencer � mettre les taxas � la position 1
	aretePos = 1;
	boot_value=0;
	
	while (string[0] == '(')   // traiter toute la chaine
	{
		a1 = 0;
		a2 = 0;
		while( string[a2] != ')')  // traiter la paire () la plus profonde
		{
			if(string[a2] == '(') a1 = a2;  // retrouver ;a parenth�se ouvrante
			a2++;
		}
		

		// a   => contient la longueur de la chaine
		// a1  => contient le debut d'un noeud � traiter
		// a2  => contient la fin d'un noeud � traiter
		// a3  => d�limite le noeud et sa longueur

		zz = a1+1;
		VertexNumber++;  // augmenter le nombre de noeuds
		boot_value=0;
		for ( ii = a1+1; ii <= a2; ii++)
		{// decortiquer cette chaine

			if (string[ii] == ':')
			{
				xx = 0;
				a3 = ii+1;

				if( string[zz] == '%')
				{ // cela veut dire que c'est un  noeud que l'on traite

					for ( jj = zz+1; jj < ii; jj++)
					{
						if(string[jj] == '|')
							break;
						string1[xx++] = string[jj]; 
					}
					temoin3=1;
					string1[xx++] = '\0';
					numero = atoi(string1);
					
					if(string[jj] == '|' ){
						boot_value=1;
						cpt_x=0;
						jj++;
						while(string[jj] != ':')
							string4[cpt_x++] = string[jj++];
						string4[cpt_x] = '\0';
					}
					
				}
				else
				{
					// on recup�re le nom du taxa

					for(jj = zz; jj < ii; jj++)
					{
						lesNoms[taxaPos-1][xx++] = string[jj];
					}
					numero = taxaPos;
					lesNoms[taxaPos-1][xx] = '\0';  // mettre la fin de chaine
					taxaPos++;  // augmenter le nombre de taxas pris
				}

			}
			else if(string[ii] == ','  || string[ii] == ')')
			{
				xx = 0;
				zz = ii +1;   // faire pointer sur le prochain noeud
				for ( jj = a3; jj < ii; jj++)
				{
					string1[xx++] = string[jj]; 
				}
				string1[xx++] = '\0';
				longueur = atof(string1);
				ARETE[aretePos++] = VertexNumber;
				ARETE[aretePos++] = numero;
				LONGUEUR[(aretePos-1)/ 2] = longueur;

				if(boot_value == 1){
					//printf("\n%d--%d : %lf (%s)",VertexNumber,numero,longueur,string4);
					boot_value=0;
				}
			}

		}

		// fin for pour traiter noeud
		//transcrire la nouvelle chaine
		xx = 0;
		for ( jj = 0; jj < (int)a1; jj++)
		{string2[xx++] = string[jj];}

		// ecrire le vertex
		//	char buffer[50];
		itoa_(VertexNumber,string1,10);
		string2[xx++] = '%';   // indiquer que c'est un noeud
		for( jj = 0; jj < (int) strlen(string1); jj++)
		{string2[xx++] = string1[jj];}

		int temoin=0;
		
		// transcrire la fin
		for( jj = a2+1; jj <= a; jj++)  // il faut voir si c'est  <= a ou c < a
		{
			if((string[jj] != ':') && (temoin==0)){
				string2[xx++] = '|';	
			}
			temoin = 1;
			//if(temoin==1)
				string2[xx++] = string[jj];
			
		}

		// remplacer string par string2 en remplacant les pointeurs
		
		tempString = string;
		string = string2;
		string2 = tempString;
		tempString = 0;
		a = xx;  // mettre la longueur � jour 

	} // fin du while pour traiter toute la string

	int root_existance = -1;
	for( jj=n;jj>0;jj--){
		strcpy(lesNoms[jj],lesNoms[jj-1]);
		if (strcmp(lesNoms[jj],"Root") == 0)
			root_existance = jj;
	}
	//printf("\nRoot_existance = %d",root_existance);
	for( jj=1;jj<=n;jj++){
		//printf("\n(%d) - %s",jj,lesNoms[jj]); 
	}
	
	
	ARETE[aretePos++] = 0;
	ARETE[aretePos++] = 0;


	for(i=1;i<=2*n-3;i++){
		LONGUEUR[i-1] = LONGUEUR[i];
		if(LONGUEUR[i-1] < 5*epsilon){
			LONGUEUR[i-1] = 5*epsilon;
		}
	}
	for(i=1;i<=2*(2*n-3);i++){
		ARETE[i-1] = ARETE[i];
	}
	
	//== recherche des branches connexes a celle de la racine;
	if( root_existance > 0){
		int noeud_interne = -1;
		for(i=1;i<=2*n-3;i++){
			if(ARETE[2*i-1] == root_existance)
				noeud_interne = ARETE[2*i-2];
			if(ARETE[2*i-2] == root_existance)
				noeud_interne = ARETE[2*i-1];
		}
		printf("\n[%d,%d]",noeud_interne,root_existance);
		for(i=1;i<=2*n-3;i++){
			if((ARETE[2*i-1] != root_existance) && (noeud_interne == ARETE[2*i-2])){
				LONGUEUR[i-1] = 50;
				printf("\n[%d,%d]",noeud_interne,ARETE[2*i-1]);
			}
			if((ARETE[2*i-2] != root_existance) && (noeud_interne == ARETE[2*i-1])){
				LONGUEUR[i-1] = 50;
				printf("\n[%d,%d]",noeud_interne,ARETE[2*i-2]);
			}
		}
		
		
	}
	
	
	//=== on teste si il y a un noeud de degre 2 et un noeud de degre 1

	//printf("\n");
	for(i=1;i<=2*n-3;i++){
		//printf("\n%d-%d --> %lf",ARETE[2*i-1],ARETE[2*i-2],LONGUEUR[i-1]);
	}
    //printf("\n");
  
	int * tableau = (int*)malloc(((2*n)+1)*sizeof(int));
	int deg2=-1,deg1=-1;
	for(i=1;i<=2*n;i++)
		tableau[i-1] = 0;
	for(i=1;i<=2*n-3;i++){
		tableau[ARETE[2*i-1]]++;
		tableau[ARETE[2*i-2]]++;
	}
//	printf("\n");
	i=n+1;
	while((tableau[i] > 0) && ((i+1)<(2*n))){
		//printf("\n%d-->%d",i,tableau[i]);
		if(tableau[i] == 2) deg2 = i;
		if(tableau[i] == 1) deg1 = i;
		i++;
	}
	//printf("\ndeg2=%d , deg1=%d",deg2,deg1);
	int pos_racine=-1;
	
	for(i=1;i<=2*n-3;i++){
		if(ARETE[2*i-1] == deg1){ ARETE[2*i-1] = deg2; pos_racine=i; break;} 
		if(ARETE[2*i-2] == deg1){ ARETE[2*i-2] = deg2; pos_racine=i; break;} 
	}
    if(pos_racine != -1){
        LONGUEUR[pos_racine-1] = 100;
    }
	/*printf("\n");
	for(i=1;i<=2*n-3;i++){
		printf("\n%d-%d --> %lf",ARETE[2*i-1],ARETE[2*i-2],LONGUEUR[i-1]);
	}
  printf("\n");
*/

	free(string);
	free(string1);
	free(string2);
	free(string4);
	free(tempString);
	free(tableau);
	//free(string3);

	(*kt) = 2*n-3 - (*kt);
    
    
	return pos_racine;

}

//===============================================================================================
//
//===============================================================================================
int lectureNewick1(const char * newick, long int * ARETE, double * LONGUEUR, char ** lesNoms)
{
	int n;                                     

	// Ce sous programme permet de lire un arbre au format newick et de le transcrire dans
	// des vecteurs arete-longueur commencant � 1
	// ATTENTION: les noms commencent � 0
	// 7 octobre 2004
	// Elmaestro

	// TODO: Add your command handler code here
	int i,j,j1,k, a, a1, a2,a3, VertexNumber,numero;
	char symbol, *string, *string1, *string2/* *string3,c ,**Et*/;
	int taxaPos; // le nombre de taxas recup�r�
	int aretePos; // le nombre d'aretes recup�r�
	char symbolOld =' ';
	int zz, xx,jj,ii;
	double longueur;
	char * tempString;
	int cpt=0;
	//long int *ARETE;
	//double *LONGUEUR; 	

	//Correctness of the Newick format verification
	i=0;
	n = 0;

	do
	{
		symbol = newick[cpt++];
		if (symbol==':') i++;
		if (symbol == ':' && symbolOld != ')') n++;
		symbolOld = symbol;
	}  while(symbol != '\0');

	cpt=0;

	if (i == 0) { printf("Incorrect Newick file format. Edge lengths must be indicated after a ':' characters."); exit(FAIL);}

	j=0;
	do{ 
	   symbol=newick[cpt++]; 
	   if (symbol=='(') j++; 
	}while(symbol != '\0');
	
	cpt=0;	
	j1=0;
	
	do{
		symbol=newick[cpt++];
		if (symbol==')') j1++; 
	}while(symbol != '\0');  
	
	cpt=0;

	// verification des arit�s de l'arbre
	if (j1!=j) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); exit(FAIL);}
	//else if (j!=n-2) { printf("Incorrect Newick file format. Only trees with vertices of degree 1 and 3 are allowed by T-REX."); fclose (data); exit(FAIL);}

	k=0;

	do{ 
		symbol=newick[cpt++];
		if (symbol==',') k++; 
	}while(symbol != '\0');   
	
	cpt=0;
	//if (k!=(n-1)) { printf("Incorrect Newick file format. Number of objects must be equal to number of commas plus 1."); fclose (data); exit(FAIL);}

	a=0;

	do{ 
		symbol=newick[cpt++];
		if (symbol==';') a++; 
	}while(symbol != '\0');    
	cpt=0;

	if (a==0) { printf("Incorrect Newick file format. Newick string must be followed by a ';' character."); exit(FAIL);}
	else if (a>1) { printf("Incorrect Newick file format. Newick string must contain (in the end) only one ';' character."); exit(FAIL);}

	a=0;
	do{
		symbol=newick[cpt++];a++;
	}while(symbol == ' '); 
	
	cpt=0;

	if (symbol!='(') { printf("Incorrect Newick file format. Newick string must begin with a '(' character."); exit(FAIL);}  

	a=0;
	
	do{ symbol=newick[cpt++];
		if (symbol=='%') a++; 
	}while(symbol != '\0');

	cpt=0;

	if (a>0) { printf("Incorrect Newick file format. Newick string cannot contain \'%%\' character."); exit(FAIL);}

	do
	{
		symbol=newick[cpt++];
		if ((symbol=='(')||(symbol==','))		
		{ 
			symbol=newick[cpt++]; a=0;
			if ((symbol!='(')&&(symbol!=',')&&(symbol!=';')&&(symbol!=':'))		   
			{ 
				cpt--; 
				do{
					symbol=newick[cpt++];a++;
				}while ((symbol!=':')&&(symbol!='\0'));  
			}
			else cpt--;
			if (a>50) { printf("Incorrect Newick file format. Names of objects must not exceed 50 characters.");  exit(FAIL);}   
		}	
	}while(symbol != '\0');
	cpt=0;

	string = (char*)malloc((1000*n) * sizeof(char));
	string2 = (char*)malloc((1000*n) * sizeof(char));
	//string3 = (char*)malloc((1000*n) * sizeof(char));
	string1 = (char*)malloc((2000) * sizeof(char));

	if ((string == NULL)||(string1 == NULL)||(string2 == NULL)/*||(string3 == NULL)*/)
	{ printf("Input data are too large or not a correct Newick file chosen "); printf("version 2 ");exit(FAIL);}

	a=0;

	do{		
		symbol=newick[cpt++];
		if ((symbol!=' ')&&(symbol!='\n')&&(symbol!='\t')) { string[a++]=symbol; } 
	}while(symbol !='\0'); 

	k=0; VertexNumber=n;
	//a1 = 0;
	//a2 = 0;
	taxaPos =1;    // nous allons commencer � mettre les taxas � la position 1
	aretePos = 1;
	while (string[0] == '(')   // traiter toute la chaine
	{
		a1 = 0;
		a2 = 0;
		while( string[a2] != ')')  // traiter la paire () la plus profonde
		{
			if(string[a2] == '(') a1 = a2;  // retrouver ;a parenth�se ouvrante
			a2++;
		}
		

		// a   => contient la longueur de la chaine
		// a1  => contient le debut d'un noeud � traiter
		// a2  => contient la fin d'un noeud � traiter
		// a3  => d�limite le noeud et sa longueur

		zz = a1+1;
		VertexNumber++;  // augmenter le nombre de noeuds
		for ( ii = a1+1; ii <= a2; ii++)
		{// decortiquer cette chaine

			if (string[ii] == ':')
			{
				xx = 0;
				a3 = ii+1;

				if( string[zz] == '%')
				{ // cela veut dire que c'est un  noeud que l'on traite

					for ( jj = zz+1; jj < ii; jj++)
					{
						string1[xx++] = string[jj]; 
					}
					string1[xx++] = '\0';
					numero = atoi(string1);
				}
				else
				{
					// on recup�re le nom du taxa

					for(jj = zz; jj < ii; jj++)
					{
						lesNoms[taxaPos-1][xx++] = string[jj];
					}
					numero = taxaPos;
					lesNoms[taxaPos-1][xx] = '\0';  // mettre la fin de chaine
					taxaPos++;  // augmenter le nombre de taxas pris
				}

			}
			else if(string[ii] == ','  || string[ii] == ')')
			{
				xx = 0;
				zz = ii +1;   // faire pointer sur le prochain noeud
				for ( jj = a3; jj < ii; jj++)
				{
					string1[xx++] = string[jj]; 
				}
				string1[xx++] = '\0';
				longueur = atof(string1);
				ARETE[aretePos++] = VertexNumber;
				ARETE[aretePos++] = numero;
				LONGUEUR[(aretePos-1)/ 2] = longueur;
			}

		}

		// fin for pour traiter noeud
		//transcrire la nouvelle chaine
		xx = 0;
		for ( jj = 0; jj < (int)a1; jj++)
		{string2[xx++] = string[jj];}

		// ecrire le vertex
		//	char buffer[50];
		itoa_(VertexNumber,string1,10);
		string2[xx++] = '%';   // indiquer que c'est un noeud
		for( jj = 0; jj < (int) strlen(string1); jj++)
		{string2[xx++] = string1[jj];}

		// transcrire la fin
		for( jj = a2+1; jj <= a; jj++)  // il faut voir si c'est  <= a ou c < a
		{string2[xx++] = string[jj];}

		// remplacer string par string2 en remplacant les pointeurs
		
		tempString = string;
		string = string2;
		string2 = tempString;
		tempString = 0;
		a = xx;  // mettre la longueur � jour 

	} // fin du while pour traiter toute la string

	for( jj=n;jj>0;jj--)
		strcpy(lesNoms[jj],lesNoms[jj-1]);

	ARETE[aretePos++] = 0;
	ARETE[aretePos++] = 0;


	for(i=1;i<=2*n-3;i++){
		LONGUEUR[i-1] = LONGUEUR[i];
	}
	for(i=1;i<=2*(2*n-3);i++){
		ARETE[i-1] = ARETE[i];
	}


	free(string);
	free(string1);
	free(string2);
	//free(string3);

	return n;

}

void filtrerMatrice(double **dissSpecies, double **dissGene, char **nomsSpecies, char **nomsGene,int nbSpecies, int nbGene, char * fichier){

		int i,j,temoin;
		/*FILE * out;
		
		if((out = fopen(fichier,"a+")) == NULL){
			printf("Can't open %s",fichier);
			exit(-1);
		}*/
		
		//fprintf(out,"\n=>Filtre");
		
		
		for(i=1;i<=nbSpecies;i++){
			//printf("\n%s -",nomsSpecies[i]);
			temoin = 0;
			for(j=1;j<=nbGene;j++){
				//printf(" %s",nomsGene[j]);
				if(strcmp(nomsSpecies[i],nomsGene[j])==0)
					temoin=1;
			}
			if(temoin == 0){
				//fprintf(out,"\nSpecies:%s",nomsSpecies[i]);
				for(j=1;j<=nbSpecies;j++){
						dissSpecies[i][j] = dissSpecies[j][i] = -1;
				}
				strcpy(nomsSpecies[i],"");
			}
		}
		//printf("==============");
		for(i=1;i<=nbGene;i++){
			//printf("\n%s -",nomsGene[i]);
			temoin = 0;
			for(j=1;j<=nbSpecies;j++){
			//	printf(" %s",nomsSpecies[j]); 	
				//if(strlen(nomsSpecies[j]) > 1) 
					if(strcmp(nomsSpecies[j],nomsGene[i])==0)
						temoin=1;
			}
			if(temoin == 0){
				//fprintf(out,"\nGene:%s",nomsGene[i]);
				for(j=1;j<=nbGene;j++){
					dissGene[i][j] = dissGene[j][i] = -1;
				}
				strcpy(nomsGene[i],"");
			}
		}

		//fclose(out);
}

//== ecrire la matrice dans le fichier output
int ecrireMatrice(double **mat,const char *outfile,int taille,char **noms){
	int i,j,finalTaille;
	FILE *out;
	
	finalTaille=0;

	/*for(i=1;i<=taille;i++)
		if(mat[1][i] != -1)
			finalTaille = finalTaille+1;*/
	for(i=1;i<=taille;i++)
		if(strcmp(noms[i],"")!=0)
			finalTaille = finalTaille+1;
	if(finalTaille < 3){
		// printf("There are less than 3 same species !");
		return -1;
	}
	if((out=fopen(outfile,"w+"))==NULL){
		printf("cannot create output file !!");
		printf("==1==>%s<===",outfile);
		exit(-1);
	}
	else{
		fprintf(out,"%d",finalTaille);
		//printf("%d",finalTaille);
		for(i=1;i<=taille;i++){
			if(strcmp(noms[i],"") != 0){//if(strlen(noms[i]) > 1){
				fprintf(out,"\n%s",noms[i]);
				//printf("\n%s",noms[i]);
				for(j=1;j<=taille;j++)
					if(mat[i][j] != -1){
						//printf(" %lf",mat[i][j]);
						fprintf(out," %lf",mat[i][j]);
					}
			}
			else{
				//printf("\nspecies : 1 cas");
			}
		}
		fclose(out);
	}
	return finalTaille;
}

//== ajouter la matrice de gene dans le fichier input
void ajouterMatriceGene(double **mat,const char *outfile,int taille,char **noms){
	int i,j;
	FILE *out;
	if((out=fopen(outfile,"a"))==NULL){
		printf("cannot create output file !!");
		printf("==2==>%s<===",outfile);
		exit(-1);
	}
	else{
		fprintf(out,"\n",taille);
		for(i=1;i<=taille;i++){
			if(strcmp(noms[i],"") != 0){  //if(strlen(noms[i]) > 1){
				fprintf(out,"\n%s",noms[i]);
				for(j=1;j<=taille;j++)
					if(mat[i][j] != -1)
						fprintf(out," %lf",mat[i][j]);
			}
			else{
				//printf("\ngene : 1 cas");
			}
		}
		fclose(out);
	}
}

/************************************************************************
* Cette procedure trie la 2eme matrice lue en fonction
* du nom des especes de la 1ere matrice
* DISS: matrice du gene
* NomsDISS: tableau des noms des especes de la matrice de gene
* NomsADD: tableau contenant le noms des especes de la matrice d'especes 
*************************************************************************/
void TrierMatrices(double **DISS,char **NomsDISS,char **NomsADD,int n)
{
	int trouve, ligne,colonne,i,j;
	double ** DISS_;
	char   ** NomsDISS_;
	char noms[50];
	int * table;
	

	table = (int *) malloc((n+1)*sizeof(int));
	DISS_ = (double **) malloc((n+1)*sizeof(double*));
	NomsDISS_ = (char **) malloc((n+1)*sizeof(char*));

	for (i=0;i<=n;i++)
	{
		DISS_[i]=(double*)malloc((n+1)*sizeof(double));
		NomsDISS_[i] = (char*) malloc((n+1)*sizeof(15));
	}

	for(ligne = 1;ligne<=n;ligne++)
	{
		strcpy(noms,NomsADD[ligne]);
		trouve = 0;
		for(colonne = 1;colonne<=n;colonne++)
		{
			if(strcmp(noms,NomsDISS[colonne])==0)
			{
				trouve = 1;
				table[ligne] = colonne;
				strcpy(NomsDISS_[ligne],noms);
			}
		}
		if(trouve==0)
		{
			printf("\n%s %s",noms,"is not in the gene matrix.This program must stop");
			exit (-1);	
		}
	}

	for(i=1;i<=n;i++)
		for(j=1;j<=i;j++) DISS_[i][j] = DISS_[j][i] = DISS[table[i]][table[j]];

	for(i=1;i<=n;i++)
	{
		strcpy(NomsDISS[i],NomsDISS_[i]);
		for(j=1;j<=n;j++) DISS[i][j] = DISS_[i][j];
	}
//exit(0);
	for (i=0;i<=n;i++)
	{
		free(DISS_[i]);
		free(NomsDISS_[i]);
	}
	free(DISS_);
    free(NomsDISS_);
	free(table);
}

/*=======================================
//
//=======================================*/
int compare(int * B, int *BI,int n){

	int somme = 0,i;
	int *tab = (int*)malloc((n+1)*sizeof(int));

	for(i=1;i<=n;i++)
		if(B[i] == BI[i])
			tab[i] = 1;
		else
			tab[i] = 0;
	for(i=1;i<=n;i++)
		somme = somme + tab[i];

	free(tab);

	return somme;

}

/*=======================================
//
//=======================================*/
void CalculerBootARETE(int *BSARETE,double **DISTG,double **DISS,int n){

	int m,mI,i,j,val;
	int **B=(int **) malloc((2*n-1)*sizeof(int*));
	int **BI=(int **) malloc((2*n-1)*sizeof(int*));
	int *PLACE=(int *) malloc((2*n-2)*sizeof(int));
	int *PLACEI=(int *) malloc((2*n-2)*sizeof(int));

	for (i=0;i<=2*n-2;i++){
		B[i]=(int *) malloc((2*n-2)*sizeof(int));
		BI[i]=(int *) malloc((2*n-2)*sizeof(int));
	}

	//TrierMatrices(DISS,NomsGene,NomsG,n);

	m=Bipartition_Table(DISTG,B,PLACE,n);
	mI=Bipartition_Table(DISS,BI,PLACEI,n);

	for(i=1;i<=2*n-3;i++){
		for(j=1;j<=2*n-3;j++){
			val = 0;
			val = compare(B[i],BI[j],n);

			if(val == n){
				BSARETE[i] = BSARETE[i] + 1;
			}	
		}
	}

	/*	printf("\n");
	for(i=1;i<=2*n-3;i++){
	printf("\n%d :",i );
	for(j=1;j<=n;j++){
	printf("%d",B[i][j]);		
	}
	}*/
	for(i=0;i<=2*n-2;i++){
		free(B[i]); free(BI[i]);
	}
	free(B);free(BI);free(PLACE);free(PLACEI);
}


void trouver_3_noeuds(double ** Adjacence,int sommet,int n,int tab_noeuds[]){
	int j=1;
	for(int i=1;i<=2*n-2;i++){
		if(Adjacence[sommet][i] < INFINI){
			tab_noeuds[j++] = i;
		}
	}
}
void trouver_fils(double **Adjacence,int sommet,int parent,int tab_fils[],int n){
	
	for(int i=1;i<=2*n-2;i++){
		if(Adjacence[sommet][i] < INFINI){
			if(i != parent){
				if(i<=n){
					tab_fils[0] += 1;
					tab_fils[tab_fils[0]] = i;	
				}
				else{
					trouver_fils(Adjacence,i,sommet,tab_fils,n);
				}
			}
		}
	}
}

void creerNoeudDT(int position,double ** DIST,struct DescTree *DT, int * tab_fils_1, int * tab_fils_2){
		
		int j;
		DT[position].Matrice = (double**)malloc((tab_fils_1[0] + tab_fils_2[0] + 1)*sizeof(double*));
		for(j=0;j<tab_fils_1[0] + tab_fils_2[0] + 1;j++){
			DT[position].Matrice[j] = (double*)malloc((tab_fils_1[0] + tab_fils_2[0] + 1)*sizeof(double));
		}
		DT[position].Tableau = (int*) malloc( (tab_fils_1[0] + tab_fils_2[0] + 1)*sizeof(int));
		DT[position].Tableau[0] = tab_fils_1[0] + tab_fils_2[0];
		for(j=1;j<=tab_fils_1[0];j++) DT[position].Tableau[j] = tab_fils_1[j];
		for(j=1;j<=tab_fils_2[0];j++) DT[position].Tableau[j+tab_fils_1[0]] = tab_fils_2[j];
		sortIntTab(DT[position].Tableau,1,DT[position].Tableau[0]);
		for(j=1;j<=DT[position].Tableau[0];j++){
			for(int k=1;k<=DT[position].Tableau[0];k++){
				DT[position].Matrice[j][k] = DIST[DT[position].Tableau[j]][DT[position].Tableau[k]];
			}
		}
		DT[position].nbSommet = DT[position].Tableau[0];
}

void deleteBipartitionSansRacine(struct DescTree * DT,int taille){
	
	int i,j;
	
	for(i=3*(taille+1);i<=3*(2*taille-2);i++){
		
		for(j=0;j<DT[i].Tableau[0];j++){
			free(DT[i].Matrice[j]);
		}
		free(DT[i].Matrice);
		free(DT[i].Tableau);
	}
	free(DT);
}


void printNoeudDT(struct DescTree *DT,int position,int n1,int n2){		
		int j,k;
		printf("\ncouple (%d,%d) position (%d): ",n1,n2,position);
		for(j=1;j<=DT[position].nbSommet;j++)
			printf("%d ",DT[position].Tableau[j]);
	/*for(j=1;j<=DT[position].nbSommet;j++){
			printf("\n%d\t",DT[position].Tableau[j]);
			for(k=1;k<=DT[position].nbSommet;k++){
				printf("%lf\t",DT[position].Matrice[j][k]);
			}
		}*/
}

//===============================================================================================================================
//
//===============================================================================================================================
void RechercherBipartitionSansRacine(long int *ARETE,double ** DIST,double ** Adjacence2,struct DescTree *DT,const int n,int kt)
{
	/*/===== variables =====*/
	int i,j,k,l,taille,plus,ns=2*n-2;
	double ** Adjacence =(double**)malloc((2*n+1)*sizeof(double*));
	int tab_noeuds[4];
	int tab_fils_1[n];
	int tab_fils_2[n];
	int tab_fils_3[n];
	int tab_fils[4][n];
	struct TNoeud * arbre;
	
	for (i=0;i<=2*n;i++)
		Adjacence[i]=(double*)malloc((2*n+1)*sizeof(double));

	
	
	//== copie de la matrice d'adjacence
	for(i=0;i<=2*n-2;i++){
		//if(i>0) printf("\n%d\t",i);
		for(j=0;j<=2*n-2;j++){
			Adjacence[i][j] = Adjacence2[i][j];
			// if(j>0 && i>0)
				// if(Adjacence2[i][j] >= INFINI){
					// printf("INF\t");
				// }
				// else{
					// printf("%1.1lf\t",Adjacence2[i][j]);
				// }
		}
	}
	
	for(i=n+1;i<=ns;i++){
		//printf("\nNoeud interne %d (nb sommets=%d)",i,ns);	
		//= trouver les 3 noeuds adjacents		
		trouver_3_noeuds(Adjacence,i,n,tab_noeuds);
		//printf("\nnoeud %d : %d,%d,%d",i,tab_noeuds[1],tab_noeuds[2],tab_noeuds[3]);	
		
		//= trouver les feuilles pour chaque fils
		tab_fils[1][0] = tab_fils[2][0] = tab_fils[3][0] = 0;
		
		for(j=1;j<=3;j++){
			if(tab_noeuds[j] <= n){ tab_fils[j][0] = 1; tab_fils[j][tab_fils[j][0]] = tab_noeuds[j];}
			trouver_fils(Adjacence,tab_noeuds[j],i,tab_fils[j],n);
		//	printf("\nfils %d : ", tab_noeuds[j]);
			//for(k=1;k<=tab_fils[j][0];k++)
				//printf("%d ",tab_fils[j][k]);
		}
			
		//== ensemble de feuilles possibles
		creerNoeudDT(3*i-2,DIST,DT,tab_fils[1],tab_fils[2]);
		//printNoeudDT(DT,3*i-2,tab_noeuds[1],tab_noeuds[2]);
		creerNoeudDT(3*i-1,DIST,DT,tab_fils[1],tab_fils[3]);
		//printNoeudDT(DT,3*i-1,tab_noeuds[1],tab_noeuds[3]);
		creerNoeudDT(3*i,DIST,DT,tab_fils[2],tab_fils[3]);
		//printNoeudDT(DT,3*i,tab_noeuds[2],tab_noeuds[3]);
			
	}	
			
			
	//== liberation de la m�moire
	for (i=0;i<=2*n;i++)
		free(Adjacence[i]);
	free(Adjacence);
}


