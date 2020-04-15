#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<sys/time.h>
#define SIZE 100

int ngauss=2;
int magmax=10;
int nmax;
int ntmax;
double dmin;
double dicmin;
int nicmax;
int nsize=8;
double coefn1[8][3]={{1.0,1.0,1.0},{-1.0,1.0,1.0},{-1.0,-1.0,1.0},{1.0,-1.0,1.0},
				{1.0,1.0,-1.0},{-1.0,1.0,-1.0},{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0}};
double coefe1[12][3]={{0.0,1.0,1.0},{0.0,-1.0,1.0},{0.0,-1.0,-1.0},{0.0,1.0,-1.0},
					 {1.0,0.0,1.0},{1.0,0.0,-1.0},{-1.0,0.0,-1.0},{-1.0,0.0,1.0},
					 {1.0,1.0,0.0},{-1.0,1.0,0.0},{-1.0,-1.0,0.0},{1.0,-1.0,0.0}};
int inode1[12][2]={{1,0},{2,3},{6,7},{5,4},{3,0},{7,4},{6,5},{2,1},{4,0},{5,1},{6,2},{7,3}};
int inode11[12][2]={{1,0},{2,3},{6,7},{5,4},{3,0},{7,4},{6,5},{2,1},{4,0},{5,1},{6,2},{7,3}};
int ir[8]={0,1,2,0,1,2,0,1};
int ixyze1[12]={0,0,0,0,1,1,1,1,2,2,2,2};

double gettime(void)
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return ((double)tv.tv_sec+(double)tv.tv_usec*1e-6);
}


int rddata(int**,int,int*,double**,int,int*,				\
		   int**,int*,int**,int*,int**,int*,int**,			\
		   double**,int*,int**,int*,double**,int*,			\
		   int**,int*,int*,double*,int*,int*,				\
		   double*,double*,int*,int*,double*,				\
		   double**,int*,double*,double*,int,int,			\
		   int**,int**,int**,int**,int,int,int,				\
		   int*,int**,double**,int*,int**,double,int,FILE*,int,int*);
int gauint(double*,double*,int,int,int);
int rdmesh(int**,int,int,int*,double**,int,int,int*,	\
		   double,int,int*,int,int,int,int,FILE*);
int edgecl(int**,int*,int,int,int**,int,int*,int*,		\
		   int**,int**,int**,int**,int,int);
int shell(int**,int,int,int);
int bubble(int**,int,int,int);
int rdpre(int**,double**,int,int*,int**,int,int*,		\
		  double**,int,int*,int**,int,int*,				\
		  int*,double*,int*,int*,double*,double*,		\
		  int*,int*,double*,int,int,int*,int**,			\
		  double**,int*,int**,double**,int,int*);
int dircal(int**,int,int*,double**,int,int*,double,		\
		   int**,int*,int**,double**,int,int*,int**,	\
		   int*,int**,double**,int*,int*,int**,int**,	\
		   int**,int,int,int*,int**,double**,int,		\
		   int*,int,int);
int initao(int*,double**,int*);
double fdiri(int,int,double*,int**,double**,int,int,int*,int*,int*);
double fjacob(int,int,int,double*,double*,double**,			\
			  double**,double**,double**,int**,double**,	\
			  int,int,int*,int*,int*);
int fem(int**,int,int*,double**,int,int*,int**,int*,		\
		int**,int*,int**,int*,int**,double**,int*,int**,		\
		int*,double**,int*,int**,int*,int*,					\
		double*,int*,int*,double*,double*,int*,int*,		\
		double*,double**,int*,double*,double*,int,			\
		int,int**,int**,int**,int**,int,int,int,int*,			\
		int**,double**,int*,int**,int*);
int ptifa(int**,int*,int**,int**,int*,int**,int*,int,
		  int*,int,int,int**,int*,int**,int*,
		  int,int*,int**,int**,int*,int,FILE*,int*);
int bandcl(void);
int iccgnb(int,int**,int*,int**,int*,int**,int**,int,
		   int*,int**,int*,int**,int**,int**,int**,int,
		   int,int,int*,int**,int,int**,int*,FILE*);
int initst(int,int*,double**);
int wtcal(double*,int*,double*,int,double*);
int setdir(int**,int*,int**,double**,int*,double**,int*,
		   int,double*,int*,int**,int*,
		   int**,double**,int,int*);
int setcur(double**,double**,int*,int*,int,int*,double*);
int gujcal(int**,int,int*,double**,int,int*,int**,
		   int*,int**,int*,int**,int*,double**,
		   double*,double*,int,int**,double**,int,int,double**,int*);
int gumcal(int**,int,int*,double**,int,int*,int**,
		   int*,int*,int**,double**,int**,int*,double**,
		   double*,double*,int,int**,double**,int,int,int**,int**,
		   int**,int*,double**,int*);
int initnr(int*,int*,double**,double**,double**);
int gss(void);
int icghfi(int**,int,int*,double**,int,int*,int**,int*,
		   int**,int*,double**,int*,double*,int**,double**,double**,int*,
		   double*,double*,int,int**,int**,int*,int**,int*,double**,double**,
		   int,double*,int*,int*,int*,FILE*);
int icgcl(int**,int,int*,double**,int,int*,int**,int*,
		  int**,int*,double**,int*,int**,double**,int*,double*,double*,
		  int,int**,int**,int*,int**,int*,double**,double**,int,
		  int,int*,int*,int*,double**,int);
int iccg(int*,int*,double**,double**,double**,int**,int**,int**,int*,double*,
		 int,int,int,FILE*,int,double**,double**,double**,double**);
int banda(int**,int,int*,double**,int,int*,int**,
		  int*,int**,double**,int*,double**,double*,int*,double**,double**,
		  int,double*,int**,int*);
int conver(int,int*,double*,double**,double*,int,int,int*,FILE*);
int eddy(int**,int,int*,double**,int,int*,int**,int*,
		 int**,int*,double**,int*,double*,double**,double**,int*,int,
		 double*,int*,FILE*);
int atoao(double**,double**,int*);
int m0cal(int**,int,int*,double**,int,int*,
		  double**,int*,int,double*,int,int,int);
int nonline(int,int,int*,double**,int**,int*,int*,double*,double*,
		    int,double*,int**,double**,int,int*,double**);
double fjacob(int,int,int,double*,double*,double**,			\
			  double**,double**,double**,int**,double**,	\
			  int,int,int*,int*,int*);
int shell(int**,int,int,int);
int bubble(int**,int,int,int);
int symsca(int,int,double**,double**,double**,double**,int**,int**,int**,int);
int inchde(int,int,double**,double**,int**,int**,int**,int,FILE*);
int forbac(int,int,double**,double**,double**,int**,int**,int**);
int forsub(int,int,double**,double**,double**,int**,int**,int**);
int bacsub(int,int,double**,double**,int**,int**,int**);

int multi(int,int,double**,double**,double**,int**,int**,int**,double);
int maxval(int,double**,double*);
int produc(int,double**,double**,double*);
int bhodf(double,double*,double*);
int bhodf21(double,double*,double*);
int isus15c(double,double*,double*);
int src(double,double*,double*);
int steel(double,double*,double*);
int ss400(double,double*,double*);
int rep_inchde(int,int,double**,double**,int**,int**,int**);
int main(void)
{
	int nelem,nvert,npoint;
	int ndim,nedges,nblock,ndiri;
	int nperi,nboun,nblok,nmate;
	int ncur,nstep,nedge,nsum;
	int ndirif,nfbun,nusef;
	int i;
	int ndebug;
	int nmag,**edge;
	char string[SIZE];
	int **nod,**idiri;
	int **iperi,**inb,**imate;
	int **icur,**imag;
	int icga;
	int **idirif,ifsin;
	int **iwh,**iwcic,**iwg,**iwgu;
	double freq,fact;
	double **xyz;
	double **amate;
	double **amag,**ao,**anb;
	double *xi,*w,dt;
	FILE *fbun;



	printf("\n");
    printf("###################################################\n");
    printf("  *     * ******* *     *    *    ******  *******  \n");
    printf("  *     * *        *   *    * *   *     *    *     \n");
    printf("  *     * *         * *    *   *  *     *    *     \n");
    printf("  ******* *****      *    *     * *     *    *     \n");
    printf("  *     * *         * *   ******* *     *    *     \n");
    printf("  *     * *        *   *  *     * *     *    *     \n");
    printf("  *     * ******* *     * *     * ******     *     \n");
    printf("###################################################\n");
    printf("                                    by iccg5       \n");
    printf("                                    by M.Shimizu   \n");
    printf("                           modified by M.Miura     \n");
    printf("                           modified by T.Imai      \n");
    printf("                           modified by H.Fusayasu  \n");
    printf("                           modified by K.Oro       \n");
    printf("                           modified by T. Torii    \n");
    printf("                           modified by H. Inoue    \n");
	printf("                         translated by T.Okimura   \n");


	nfbun=10;
	double t0=gettime();
	if(nfbun>0)
	{
		fbun=fopen("bun.hex","r");
		
		fgets(string,SIZE,fbun);
		fgets(string,SIZE,fbun);
		fscanf(fbun,"%d",&npoint);
		fgets(string,SIZE,fbun);
		fgets(string,SIZE,fbun);
		fscanf(fbun,"%d",&nelem);
		fgets(string,SIZE,fbun);
		fgets(string,SIZE,fbun);
		fscanf(fbun,"%d",&nblock);
		fgets(string,SIZE,fbun);
		fgets(string,SIZE,fbun);
		fscanf(fbun,"%le",&fact);


		printf("npoint=%d\n",npoint);
		printf("nelem=%d\n",nelem);
		printf("nblock=%d\n",nblock);
		printf("fact=%le\n",fact);
	}
	ntmax=npoint*4;
	nsum=npoint*4*nsize/4;
	printf("nsum=%d\n",nsum);
	printf("malloc start!\n");

	nod=(int **)malloc(sizeof(int*)*27);
	for(i=0;i<27;i++)
	{
		nod[i]=(int *)malloc(sizeof(int)*nelem);
		if(nod[i]==NULL)
		{
			goto malerrer;
		}
	}
	xyz=(double **)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		xyz[i]=(double *)malloc(sizeof(double)*npoint);
		if(xyz[i]==NULL)
		{
			goto malerrer;
		}
	}
	edge=(int **)malloc(sizeof(int*)*nelem);
	for(i=0;i<nelem;i++)
	{
		edge[i]=(int *)malloc(sizeof(int)*12);
		if(edge[i]==NULL)
		{
			goto malerrer;
		}
	}
	idiri=(int **)malloc(sizeof(int*)*25);
	for(i=0;i<25;i++)
	{
		idiri[i]=(int *)malloc(sizeof(int)*3);
		if(idiri[i]==NULL)
		{
			goto malerrer;
		}
	}
	iperi=(int **)malloc(sizeof(int*)*2);
	for(i=0;i<2;i++)
	{
		iperi[i]=(int *)malloc(sizeof(int)*npoint*4);
		if(iperi[i]==NULL)
		{
			goto malerrer;
		}
	}
	
	inb=(int **)malloc(sizeof(int*)*(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		inb[i]=(int *)malloc(sizeof(int)*100);
		if(inb[i]==NULL)
		{
			goto malerrer;
		}
	}
	anb=(double **)malloc(sizeof(double*)*(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		anb[i]=(double *)malloc(sizeof(double)*100);
		if(anb[i]==NULL)
		{
			goto malerrer;
		}
	}
	imate=(int **)malloc(sizeof(int*)*25);
	for(i=0;i<25;i++)
	{
		imate[i]=(int *)malloc(sizeof(int)*4);
		if(imate[i]==NULL)
		{
			goto malerrer;
		}
	}
	amate=(double **)malloc(sizeof(double*)*25);
	for(i=0;i<25;i++)
	{
		amate[i]=(double *)malloc(sizeof(double)*4);
		if(amate[i]==NULL)
		{
			goto malerrer;
		}
	}
	icur=(int **)malloc(sizeof(int*)*25);
	for(i=0;i<25;i++)
	{
		icur[i]=(int *)malloc(sizeof(int)*4);
		if(icur[i]==NULL)
		{
			goto malerrer;
		}
	}
	ao=(double **)malloc(sizeof(double*)*(npoint/25+1));
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		ao[i]=(double *)malloc(sizeof(double)*100);
		if(ao[i]==NULL)
		{
			goto malerrer;
		}
	}
	xi=(double *)malloc(sizeof(double)*10);
	w=(double *)malloc(sizeof(double)*10);
	iwg=(int **)malloc(sizeof(int*)*((int)(npoint/25+1)));
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		iwg[i]=(int *)malloc(sizeof(int)*100);
		if(iwg[i]==NULL)
		{
			goto malerrer;
		}
	}
	iwgu=(int **)malloc(sizeof(int*)*((int)(npoint/25+1)));
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		iwgu[i]=(int *)malloc(sizeof(int)*100);
		if(iwgu[i]==NULL)
		{
			goto malerrer;
		}
	}
	iwh=(int **)malloc(sizeof(int*)*nelem);
	for(i=0;i<nelem;i++)
	{
		iwh[i]=(int *)malloc(sizeof(int)*100);
		if(iwh[i]==NULL)
		{
			goto malerrer;
		}
	}
	iwcic=(int **)calloc(nelem,sizeof(int*));
	for(i=0;i<nelem;i++)
	{
		iwcic[i]=(int *)calloc(100,sizeof(int));
		if(iwcic[i]==NULL)
		{
			goto malerrer;
		}
	}
	imag=(int **)malloc(sizeof(int*)*10);
	for(i=0;i<10;i++)
	{
		imag[i]=(int *)malloc(sizeof(int)*2);
	}
	amag=(double **)malloc(sizeof(double*)*10);
	for(i=0;i<10;i++)
	{
		amag[i]=(double *)malloc(sizeof(double)*3);
	}
	idirif=(int **)malloc(sizeof(int*)*10);
	for(i=0;i<10;i++)
	{
		idirif[i]=(int *)malloc(sizeof(int)*3);
	}

	printf("malloc end!\n");


	ndebug=0;    /*0�ł͂Ȃ��Ƃ��f�o�b�O���[�h*/

	/*�f�[�^�̓ǂݍ���*/

	rddata( nod   ,nelem  ,&nvert ,xyz   ,npoint ,&ndim  ,	\
    edge  ,&nedges,idiri ,&ndiri ,iperi ,&nperi ,inb   ,	\
	anb   ,&nboun ,imate ,&nblok ,amate ,&nmate ,icur  ,	\
	&ncur  ,&nmax  ,&dmin  ,&icga  ,&nicmax,&dicmin,		\
	&freq  ,&ifsin ,&nstep ,&dt    ,ao    ,&nedge ,xi    ,	\
	w     ,ngauss ,ntmax  ,iwh   ,iwcic ,iwg   ,iwgu  ,		\
	nsum   ,ndebug ,magmax ,&nmag  ,imag  ,amag  ,&ndirif,	\
	idirif,fact   ,nblock,fbun,nfbun,&nusef);
	if(0)
	{
malerrer:;
		printf("malloc errer!!!!\n program ended.");
	}

	/*FEM*/
	double t1=gettime();
	fem(nod   ,nelem  ,&nvert ,xyz   ,npoint ,				\
	    &ndim  ,edge  ,&nedges,idiri ,&ndiri ,				\
		iperi ,&nperi ,inb   ,anb   ,&nboun ,				\
        imate ,&nblok ,amate ,&nmate ,icur  ,				\
		&ncur  ,&nmax  ,&dmin  ,&icga  ,					\
		&nicmax,&dicmin,&freq  ,&ifsin ,&nstep ,			\
		&dt    ,ao    ,&nedge ,xi    ,w     ,				\
		ngauss ,ntmax  ,iwh   ,iwcic ,iwg   ,				\
		iwgu  ,nsum   ,ndebug ,magmax ,&nmag  ,				\
		imag  ,amag  ,&ndirif,idirif,&nusef   );
	double t2=gettime();
	printf("\n All Computation time: tc= %9.7e[s]",(t2-t0));
	printf("\n RDD Computation time: tr= %9.7e[s]",(t1-t0));
	printf("\n FEM Computation time: tf= %9.7e[s]",(t2-t1));
	printf("FEM ending!\n");

	return 0;
}
int rddata( int **nod   ,int nelem  ,int *nvert ,double **xyz   ,int npoint ,int *ndim  ,						
    int **edge  ,int *nedges,int **idiri ,int *ndiri ,int **iperi ,int *nperi ,int **inb   ,				
	double **anb   ,int *nboun ,int **imate ,int *nblok ,double **amate ,int *nmate ,int **icur  ,		
	int *ncur  ,int *nmax  ,double *dmin  ,int *icga  ,int *nicmax,double *dicmin,		
	double *freq  ,int *ifsin ,int *nstep ,double *dt    ,double **ao    ,int *nedge ,double *xi    ,	
	double *w     ,int ngauss ,int ntmax  ,int **iwh   ,int **iwcic ,int **iwg   ,int **iwgu  ,				
	int nsum   ,int ndebug ,int magmax ,int *nmag  ,int **imag  ,double **amag  ,int *ndirif,			
	int **idirif ,double fact   ,int nblock,FILE* fbun,int nfbun,int *nusef)
{
	int nvmax=27;
	int ndmax=3;
	int negmax=36;
	int nblmax=25;
	int nbkmax=nblmax;
	int nrmax=10;
	int nmamax=nblmax;
	int ncumax=nblmax;
	int nrfmax=10;
	int ngamax=10;
	int nbfmax=npoint;
	double **adirif,**adiri;
	int i,nbounf;
	FILE *fa;

	printf("\n\n program readdata start...\n\n");

	gauint(xi    ,w     ,ngamax,ngauss,ndebug);

	printf("subroutine gaiunt ok!\n");
	 
	

    rdmesh(nod,nvmax,nelem,nvert,				\
    xyz   ,ndmax ,npoint,ndim  ,fact,			\
    negmax,nedges,								\
    nblmax,nblock,								\
    nfbun ,ndebug,fbun);

	printf("subroutine rdmesh ok!\n");
	
	adirif=(double **)malloc(sizeof(double*)*nrfmax);
	for(i=0;i<nrfmax;i++)
	{
		adirif[i]=(double *)malloc(sizeof(double)*2);
	}
	adiri=(double **)malloc(sizeof(double*)*nrmax);
	for(i=0;i<nrmax;i++)
	{
		adiri[i]=(double *)malloc(sizeof(double)*4);
	}

    rdpre(idiri, adiri, nrmax, ndiri,				\
    imate, nbkmax, nblok, amate, nmamax, nmate,		\
    icur, ncumax, ncur,								\
    nmax, dmin, icga, nicmax, dicmin,				\
    freq, ifsin,nstep, dt, ndebug,					\
    magmax,nmag, imag, amag,						\
    ndirif,idirif,adirif,nrfmax,nusef);

	printf("subroutine rdpre ok!\n");

    edgecl(nod   ,nvert ,nelem ,npoint,			\
    edge  ,negmax,nedges,nedge,					\
    iwh   ,iwcic ,iwg   ,iwgu  ,nsum  ,ndebug);

	printf("subroutine edgecl ok!\n");

    dircal(nod, nelem, nvert, xyz, npoint,		\
	ndim, fact,edge, nedges, idiri, adiri,		\
	nrmax, ndiri,iperi, nperi, inb,				\
	anb, nboun,nedge,							\
    iwh, iwcic, iwg, nsum, ndebug,				\
    ndirif,idirif,adirif,nrfmax,				\
    &nbounf,nbfmax,negmax);

	printf("subroutine dircal ok!\n");

	for(i=0;i<nrmax;i++)
	{
		free(adiri[i]);
	}
	free(adiri);

	for(i=0;i<nrfmax;i++)
	{
		free(adirif[i]);
	}
	free(adirif);

    initao(ifsin, ao, nedge);

	printf("subroutine initao ok!\n");
	fa=fopen("inb.dat","w");
	for(i=0;i<*nboun;i++)
	{
		fprintf(fa,"%d\t%e\n",inb[(int)i/100][(int)fmod(i,100)],anb[(int)i/100][(int)fmod(i,100)]);
	}
	fclose(fa);

	printf("\n\n program readdata ok!\n\n");

    return 0;
}

int gauint(double *xi,double *w,int ngamax,int ngauss,int ndebug)
{
	int i;
	if(ngauss>ngamax)
	{
		printf("ngauss is improper value.\n");
		
	}
	else if (ngauss==1)
	{
        xi[0] =   0.000000000000000;
        w[0]  =   2.000000000000000;
	}
    else if (ngauss==2)
	{
        xi[0] = - 0.577350269189626;
        xi[1] = -xi[0];
        w[0]  =   1.000000000000000;
        w[1]  =   w[0];
	}
    else if (ngauss==3)
	{
        xi[0] = - 0.774596669241483;
        xi[1] =   0.000000000000000;
        xi[2] = - xi[0];
        w[0]  =   0.555555555555556;
        w[1]  =   0.888888888888889;
        w[2]  =   w[0];
	}
    else if (ngauss==4)
	{
        xi[0] = - 0.861136311594053;
        xi[1] = - 0.339981043584856;
        xi[2] = - xi[1];
        xi[3] = - xi[0];
        w[0]  =   0.347854845137454;
        w[1]  =   0.652145154862546;
        w[2]  =   w[1];
		w[3]  =   w[0];
	}
    else if (ngauss==5)
	{
        xi[0] = - 0.906179845938664;
        xi[1] = - 0.538469310105683;
        xi[2] =   0.000000000000000;
        xi[3] = - xi[1];
        xi[4] = - xi[0];
        w[0]  =   0.236926885056189;
        w[1]  =   0.478628670499366;
        w[2]  =   0.568888888888889;
        w[3]  =   w[1];
        w[4]  =   w[0];
	}
    else if (ngauss==6)
	{
        xi[0] = - 0.932469514203152;
        xi[1] = - 0.661209386466265;
        xi[2] = - 0.238619186083197;
        xi[3] = - xi[2];
        xi[4] = - xi[1];
        xi[5] = - xi[0];
        w[0]  =   0.171324492379170;
        w[1]  =   0.360761573048139;
        w[2]  =   0.467913934572691;
        w[3]  =   w[2];
        w[4]  =   w[1];
        w[5]  =   w[0];
	}
    else if (ngauss==7)
	{
        xi[0] = - 0.949107912342759;
        xi[1] = - 0.741531185599394;
        xi[2] = - 0.405845151377397;
        xi[3] =   0.000000000000000;
        xi[4] = - xi[2];
        xi[5] = - xi[1];
        xi[6] = - xi[0];
        w[0]  =   0.129484966168870;
        w[1]  =   0.279705391489277;
        w[2]  =   0.381830050505119;
        w[3]  =   0.417959183673469;
        w[4]  =   w[2];
        w[5]  =   w[1];
		w[6]  =   w[0];
	}
    else if (ngauss==8)
	{
        xi[0] = - 0.960289856497536;
        xi[1] = - 0.796666477413627;
        xi[2] = - 0.525532409916329;
        xi[3] = - 0.183434642495650;
        xi[4] = - xi[3];
        xi[5] = - xi[2];
        xi[6] = - xi[1];
        xi[7] = - xi[0];
        w[0]  =   0.101228536290376;
        w[1]  =   0.222381034453374;
        w[2]  =   0.313706645877887;
        w[3]  =   0.362683783378362;
        w[4]  =   w[3];
        w[5]  =   w[2];
        w[6]  =   w[1];
        w[7]  =   w[0];
	}
    else if (ngauss==9)
	{
        xi[0] = - 0.968160239507626;
        xi[1] = - 0.836031107326636;
        xi[2] = - 0.613371432700590;
        xi[3] = - 0.324253423403809;
        xi[4] =   0.000000000000000;
        xi[5] = - xi[3];
        xi[6] = - xi[2];
        xi[7] = - xi[1];
        xi[8] = - xi[0];
        w[0]  =   0.081274388361574;
        w[1]  =   0.180648160694857;
        w[2]  =   0.260610696402935;
        w[3]  =   0.312347077040003;
        w[4]  =   0.330239355001260;
        w[5]  =   w[3];
        w[6]  =   w[2];
        w[7]  =   w[1];
        w[8]  =   w[0];
	}
    else if (ngauss==10)
	{
        xi[0] = - 0.973906528517172;
        xi[1] = - 0.865063366688985;
        xi[2] = - 0.679409568299024;
        xi[3] = - 0.433395394129247;
        xi[4] = - 0.148874338981631;
        xi[5] = - xi[4];
        xi[6] = - xi[3];
        xi[7] = - xi[2];
        xi[8] = - xi[1];
        xi[9] = - xi[0];
        w[0]  =   0.066671344308688;
        w[1]  =   0.149451349150581;
        w[2]  =   0.219086362515982;
        w[3]  =   0.269266719309996;
        w[4]  =   0.295524224714753;
        w[5]  =   w[4];
        w[6]  =   w[3];
        w[7]  =   w[2];
        w[8]  =   w[1];
        w[9]  =   w[0];
	}
	else
	{
		printf("Please input positive number in ngauss.\n");
	}

	if(ndebug>0)
	{
		printf("ngauss = %d\n",ngauss);
		for(i=0;i<ngamax;i++)
		{
			printf("xi(%d) = %le\n",i,xi[i]);
		}
		for(i=0;i<ngamax;i++)
		{
			printf("w(%d)  = %le\n",i,w[i]);
		}
	}
	return 0;
}

int rdmesh(int **nod,int nvmax,int nelem,int *nvert,double **xyz,int ndmax,int npoint,int *ndim,	\
	double fact,int negmax,int *nedges,int nblmax,int nblock,int nfbun,int ndebug,FILE* fbun)
{
	int i,j;
	int nd;


	*nvert=8;
	*ndim=3;
	*nedges=12;


	printf("\n\n Now is reading mesh data! \n\n");

	nd=*ndim;

	for(i=0;i<nelem;i++)
	{
/*		fscanf(fbun,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",&nod[0][i],&nod[1][i],&nod[2][i],&nod[3][i],&nod[4][i],&nod[5][i],&nod[6][i],&nod[7][i]);*/
fscanf(fbun,"%d,%d,%d,%d,%d,%d,%d,%d,",&nod[0][i],&nod[1][i],&nod[2][i],&nod[3][i],&nod[4][i],&nod[5][i],&nod[6][i],&nod[7][i]);
		if(ndebug>0)
		{
			printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",nod[0][i],nod[1][i],nod[2][i],nod[3][i],nod[4][i],nod[5][i],nod[6][i],nod[7][i]);
		}
	}
	for(i=0;i<npoint;i++)
	{
/*		fscanf(fbun,"%le\t%le\t%le",&xyz[0][i],&xyz[1][i],&xyz[2][i]);*/
		fscanf(fbun,"%le,%le,%le,",&xyz[0][i],&xyz[1][i],&xyz[2][i]);
		if(ndebug>0)
		{
			printf("%e\t%e\t%e\n",xyz[0][i],xyz[1][i],xyz[2][i]);
		}
	}
	printf("xyz reading finish!\n");

	if(nfbun>0)
	{
		fclose(fbun);
	}
	
	for(i=0;i<npoint;i++)
	{
		for(j=0;j<nd;j++)
		{
			xyz[j][i]=xyz[j][i]/fact;
		}
	}

	return 0;
}

int edgecl(int **nod,int *nvert,int nelem,int npoint,int **edge,int negmax,int *nedges,int *nedge,		\
		   int **iwork,int **iwnod,int **iwg,int **iwkp,int nsum,int ndebug)
{
	int i,j,k,nod1,nod2,nend;
	int nstart,nd,nsmall,nlarge;
	int na,nb,isum,is,isa,l;
	FILE *fa;
	fa=fopen("edge.dat","w");

	for(i=0;i<=npoint;i++)
	{
		iwg[(int)i/100][(int)fmod(i,100)]=0;
	}
	isum=0;
	for(i=0;i<nelem;i++)
	{
		for(j=0;j<*nedges;j++)
		{
			nod1=nod[inode1[j][0]][i];
			nod2=nod[inode1[j][1]][i];
			
			if(nod1<nod2)
			{
				iwork[(int)isum/100][(int)fmod(isum,100)]=nod1;
				iwork[(int)(isum+1)/100][(int)fmod(isum+1,100)]=nod2;
				iwg[(int)(nod1)/100][(int)fmod(nod1,100)]=iwg[(int)(nod1)/100][(int)fmod(nod1,100)]+1;
			}
			else
			{
				iwork[(int)isum/100][(int)fmod(isum,100)]=nod2;
				iwork[(int)(isum+1)/100][(int)fmod(isum+1,100)]=nod1;
				iwg[(int)(nod2)/100][(int)fmod(nod2,100)]=iwg[(int)(nod2)/100][(int)fmod(nod2,100)]+1;
			}
			isum=isum+2;
		}
	}

	iwkp[0][0]=0;
	for(i=1;i<=npoint;i++)
	{
		iwkp[(int)i/100][(int)fmod(i,100)]=iwkp[(int)(i-1)/100][(int)fmod(i-1,100)]+iwg[(int)(i-1)/100][(int)fmod(i-1,100)];
	}
	for(i=0;i<isum;)
	{
		isa=iwork[(int)i/100][(int)fmod(i,100)];
		is=iwkp[(int)isa/100][(int)fmod(isa,100)];
		iwnod[(int)is/100][(int)fmod(is,100)]=iwork[(int)(i+1)/100][(int)fmod(i+1,100)];
		iwkp[(int)isa/100][(int)fmod(isa,100)]=iwkp[(int)isa/100][(int)fmod(isa,100)]+1;
		i=i+2;
	}
	iwkp[0][0]=0;
	for(i=1;i<=npoint;i++)
	{
		iwkp[(int)i/100][(int)fmod(i,100)]=iwkp[(int)(i-1)/100][(int)fmod(i-1,100)]+iwg[(int)(i-1)/100][(int)fmod(i-1,100)];
	}
	*nedge=0;
	nd=*nedge;
	nend=0;
	for(i=1;i<=npoint;i++)
	{
		nstart=nend;
		nend=nstart+iwg[(int)i/100][(int)fmod(i,100)];

		bubble(iwnod,nsum,nstart,nend-1);
		nd=nd+1;
		l=1;
		iwnod[(int)(nd-1)/100][(int)fmod(nd-1,100)]=iwnod[(int)(nstart)/100][(int)fmod(nstart,100)];
		if(iwg[(int)i/100][(int)fmod(i,100)]==0)
		{
			l=0;
			nd=nd-1;
			goto edge0;
		}
		for(j=nstart+1;j<nend;j++)
		{
			if(iwnod[(int)j/100][(int)fmod(j,100)]>iwnod[(int)(nd-1)/100][(int)fmod(nd-1,100)])
			{
				nd=nd+1;
				iwnod[(int)(nd-1)/100][(int)fmod(nd-1,100)]=iwnod[(int)j/100][(int)fmod(j,100)];
				l=l+1;
			}
		}
edge0:;
		iwg[(int)i/100][(int)fmod(i,100)]=l;
	}
	*nedge=nd;

	iwkp[0][0]=1;
	for(i=1;i<npoint+1;i++)
	{
		iwkp[(int)i/100][(int)fmod(i,100)]=iwkp[(int)(i-1)/100][(int)fmod(i-1,100)]+iwg[(int)(i-1)/100][(int)fmod(i-1,100)];
	}
	for(i=0;i<nelem;i++)
	{
		for(j=0;j<*nedges;j++)
		{
			nod1=nod[inode1[j][0]][i];
			nod2=nod[inode1[j][1]][i];
			if(nod1<=nod2)
			{
				nsmall=nod1;
				nlarge=nod2;
			}
			else
			{
				nsmall=nod2;
				nlarge=nod1;
			}
			na=iwkp[(int)(nsmall)/100][(int)fmod(nsmall,100)];
			nb=iwkp[(int)(nsmall)/100][(int)fmod(nsmall,100)]+iwg[(int)(nsmall)/100][(int)fmod(nsmall,100)];
			for(k=na-1;k<nb;k++)
			{
				if(iwnod[(int)k/100][(int)fmod(k,100)]==nlarge)
				{
					goto okey;
				}
			}
			printf("\n nsmall=%d\n",nsmall);
			printf(" nlarge=%d\n",nlarge);
			for(k=na;k<nb;k++)
			{
				printf("%d\t%d\n",iwnod[(int)k/100][(int)fmod(k,100)],nlarge);
			}
			printf("---errer---can't search nod!  \n");
			return 1;
okey:;
			edge[i][j]=k;
			fprintf(fa,"%d\t%d\t%d\n",i,j,edge[i][j]);
		}
	}
	fclose(fa);

	return 0;
}

int shell(int **idata,int isize,int istart,int iend)
{
	int i,j,nins,idelta,insert;

	if(istart>=iend)
	{
		return 0;
	}
	idelta=iend-istart;
reshell:;
	idelta=idelta/2;
	for(i=istart+idelta;i<iend;i++)
	{
		insert=i;
		nins=idata[(int)insert/100][(int)fmod(insert,100)];
		j=i-idelta;
conshell:;
		if(idata[(int)j/100][(int)fmod(j,100)]<nins)
		{
			goto shellok;
		}
		insert=j;
		idata[(int)(j+idelta)/100][(int)fmod(j+idelta,100)]=idata[(int)j/100][(int)fmod(j,100)];
		j=j-idelta;
		if(j>=istart)
		{
			goto conshell;
		}
shellok:;
		idata[(int)insert/100][(int)fmod(insert,100)]=nins;
	}
	if(idelta>1)
	{
		goto reshell;
	}

	return 0;
}

int bubble(int **idata,int isize,int istart,int iend)
{
	int ibase,jmin,nmin,isrch;
	if(istart>=iend)
	{
		return 0;
	}
	for(ibase=istart;ibase<=iend;ibase++)
	{
		jmin=ibase;
		nmin=idata[(int)jmin/100][(int)fmod(jmin,100)];
		for(isrch=ibase+1;isrch<=iend;isrch++)
		{
			if(idata[(int)isrch/100][(int)fmod(isrch,100)]<nmin)
			{
				jmin=isrch;
				nmin=idata[(int)jmin/100][(int)fmod(jmin,100)];
			}
		}
		idata[(int)jmin/100][(int)fmod(jmin,100)]=idata[(int)ibase/100][(int)fmod(ibase,100)];
		idata[(int)ibase/100][(int)fmod(ibase,100)]=nmin;
	}

	return 0;
}

int rdpre(int **idiri,double **adiri,int nrmax,int *ndiri,int **imate,int nbkmax,int *nblok,		\
		  double **amate,int nmamax,int *nmate,int **icur,int ncumax,int *ncur,		\
		  int *nmax,double *dmin,int *icga,int *nicmax,double *dicmin,double *freq,				\
		  int *ifsin,int *nstep,double *dt,int ndebug,int magmax,int *nmag,int **imag,				\
		  double **amag,int *ndirif,int **idirif,double **adirif,int nrfmax,int *nusef)
{
	int i;
	char string[20];
	FILE *fpre;

	fpre=fopen("pre.dat","r");
	
	fscanf(fpre,"%d",ndiri);
	printf("%d\n",*ndiri);

	for(i=0;i<*ndiri;i++)
	{
		idiri[i][0]=0;
		idiri[i][1]=0;
	}
	for(i=0;i<*ndiri;i++)
	{
		fscanf(fpre,"%d\t%le\t%le\t%le\t%le",&idiri[i][2],&adiri[i][0],&adiri[i][1],&adiri[i][2],&adiri[i][3]);
		printf("%d\t%e\t%e\t%e\t%e\n",idiri[i][2],adiri[i][0],adiri[i][1],adiri[i][2],adiri[i][3]);

	}
	fscanf(fpre,"%d",nusef);
	if(*nusef==1)
	{
		fscanf(fpre,"%d",ndirif);
	}

	fscanf(fpre,"%d\t%d",nblok,nmate);
	printf("%d\t%d\n",*nblok,*nmate);

	for(i=0;i<*nblok;i++)
	{
		fscanf(fpre,"%d\t%d\t%d\t%d",&imate[i][0],&imate[i][1],&imate[i][2],&imate[i][3]);
		printf("%d\t%d\t%d\t%d\n",imate[i][0],imate[i][1],imate[i][2],imate[i][3]);

	}
	for(i=0;i<*nmate;i++)
	{
		fscanf(fpre,"%le\t%le\t%le\t%le",&amate[i][0],&amate[i][1],&amate[i][2],&amate[i][3]);
		printf("%e\t%e\t%e\t%e\n",amate[i][0],amate[i][1],amate[i][2],amate[i][3]);
	}

	fscanf(fpre,"%d\t%le\t%d\t%d\t%le",nmax,dmin,icga,nicmax,dicmin);
	printf("%d\t%e\t%d\t%d\t%e\n",*nmax,*dmin,*icga,*nicmax,*dicmin);

	fscanf(fpre,"%d\t%d\t%le\t%d\t%d\t%le",ncur,nmag,freq,ifsin,nstep,dt);
	printf("%d\t%d\t%e\t%d\t%d\t%e\n",*ncur,*nmag,*freq,*ifsin,*nstep,*dt);
	
	if(ncur>0)
	{
		for(i=0;i<*ncur;i++)
		{
			fscanf(fpre,"%d\t%d",&icur[i][0],&icur[i][1]);
			printf("%d\t%d\n",icur[i][0],icur[i][1]);
		}
	}
	if(nmag>0)
	{
		for(i=0;i<*nmag;i++)
		{
			fscanf(fpre,"%d\t%d\t%le\t%le\t%le",&imag[i][0],&imag[i][1],&amag[i][0],&amag[i][1],&amag[i][2]);
			printf("%d\t%d\t%e\t%e\t%e\n",imag[i][0],imag[i][1],amag[i][0],amag[i][1],amag[i][2]);
		}
	}

	fclose(fpre);

	return 0;
}

int dircal(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,double fact,			\
		   int **edge,int *nedges,int **idiri,double **adiri,int nrmax,int *ndiri,int **iperi,		\
		   int *nperi,int **inb,double **anb,int *nboun,int *nedge,int **iwelem,int **iwedge,		\
		   int **iwperi,int nsum,int ndebug,int *ndirif,int **idirif,double **adirif,int nrfmax,	\
		   int *nbounf,int nbfmax,int negmax)
{
	int nod1,nod2,nnod1,nnod2;
	int i,j,k,l,m,n,na,nb,ed,nd,jk,ja,ne;
	int ke,le,np;
	int inod,jnod,iinod,jjnod;
	int knod,lnod,kknod,llnod;
	int nwperi;
	double dx,dy,dz,h,dedge,dxy;
	double *adir;
	double zh,xmin,xmax,xp,yp,zp;
	double ymin,ymax,zmin,zmax;
	double xip,yip,zip,xjp,yjp,zjp;
	double rpi,rpj,zahyo;
	
	nnod1=nod[inode1[0][0]][0];
	nnod2=nod[inode1[0][1]][0];
	if(nnod1<nnod2)
	{
		nod1=nnod1;
		nod2=nnod2;
	}
	else
	{
		nod1=nnod2;
		nod2=nnod1;
	}
	dx=xyz[0][nod1-1]-xyz[0][nod2-1];
	dy=xyz[1][nod1-1]-xyz[1][nod2-1];
	dz=xyz[2][nod1-1]-xyz[2][nod2-1];
	h=sqrt(dx*dx+dy*dy+dz*dz);
	for(i=0;i<nelem;i++)
	{
		for(j=0;j<*nedges;j++)
		{
			nnod1=nod[inode1[j][0]][i];
			nnod2=nod[inode1[j][1]][i];

			if(nnod1<nnod2)
			{
				nod1=nnod1;
				nod2=nnod2;
			}
			else
			{
				nod1=nnod2;
				nod2=nnod1;
			}

			dx=xyz[0][nod1-1]-xyz[0][nod2-1];
			dy=xyz[1][nod1-1]-xyz[1][nod2-1];
			dz=xyz[2][nod1-1]-xyz[2][nod2-1];
			dedge=sqrt(dx*dx+dy*dy+dz*dz);

			if(dedge<h)
			{
				h=dedge;
			}
		}
	}
	h=h*0.001;

	printf("h=%e\n",h);

	for(i=0;i<nelem;i++)
	{
		for(j=0;j<*nedges;j++)
		{
			ed=edge[i][j];
			iwelem[(int)ed/100][(int)fmod(ed,100)]=i;
			iwedge[(int)ed/100][(int)fmod(ed,100)]=j;
		}
	}

	nb=0;
	*nperi=0;
	np=*nperi;
	adir=(double *)malloc(sizeof(double)*3);

	nd=*ndiri;
	for(k=0;k<*ndiri;k++)
	{
		idiri[k][0]=nb;
		jk=idiri[k][2];
		if((jk>=1)&&(jk<=3))
		{
			zh=adiri[k][0];
			xmin=zh/fact-h;
			xmax=zh/fact+h;
			adir[0]=adiri[k][1];
			adir[1]=adiri[k][2];
			adir[2]=adiri[k][3];

			for(l=0;l<*nedge;l++)
			{
				i=iwelem[(int)l/100][(int)fmod(l,100)];
				j=iwedge[(int)l/100][(int)fmod(l,100)];

				iinod=nod[inode1[j][0]][i];
				jjnod=nod[inode1[j][1]][i];
				if(iinod<jjnod)
				{
					inod=iinod;
					jnod=jjnod;
				}
				else
				{
					inod=jjnod;
					jnod=iinod;
				}

				xp=(xyz[jk-1][inod-1]+xyz[jk-1][jnod-1])/2.0;
				if((xp>=xmin)&&(xp<=xmax))
				{
					inb[(int)nb/100][(int)fmod(nb,100)]=l;
					anb[(int)nb/100][(int)fmod(nb,100)]=fdiri(i,j,adir,nod,xyz,npoint,nelem,nvert,ndim,nedges);

					if(iinod>=jjnod)
					{
						anb[(int)nb/100][(int)fmod(nb,100)]=-anb[(int)nb/100][(int)fmod(nb,100)];
					}
					nb=nb+1;
				}
			}
		}
		if((jk>=4)&&(jk<=6))
		{
			xmin = adiri[k][1] / fact - h;
			ymin = adiri[k][2] / fact - h;
			zmin = adiri[k][3] / fact - h;
			xmax = adiri[k][1] / fact + h;
			ymax = adiri[k][2] / fact + h;
			zmax = adiri[k][3] / fact + h;

			nwperi=0;
			ne=*nedge;
			for(n=0;n<ne-1;n++)
			{
				i = iwelem[(int)n/100][(int)fmod(n,100)];
				j = iwedge[(int)n/100][(int)fmod(n,100)];
            
				iinod = nod[inode1[j][0]][i];
				jjnod = nod[inode1[j][1]][i];
				if (iinod<jjnod)
				{
					inod = iinod;
					jnod = jjnod;
				}
				else
				{
					inod = jjnod;
					jnod = iinod;
				}
				xp = (xyz[0][inod-1] + xyz[0][jnod-1]) / 2.0;
				yp = (xyz[1][inod-1] + xyz[1][jnod-1]) / 2.0;
				zp = (xyz[2][inod-1] + xyz[2][jnod-1]) / 2.0;
				if (((xp>xmin)&&(xp <xmax))||((yp>ymin)&&(yp<ymax))||((zp>zmin)&&(zp<zmax)))
				{
					iwperi[(int)nwperi/100][(int)fmod(nwperi,100)] = n;
					nwperi = nwperi + 1;
				}
			}

			for(n=0;n<nwperi;n++)
			{
				m=iwperi[(int)n/100][(int)fmod(n,100)];
				i=iwelem[(int)m/100][(int)fmod(m,100)];
				j=iwedge[(int)m/100][(int)fmod(m,100)];
				iinod = nod[inode1[j][0]][i];
				jjnod = nod[inode1[j][1]][i];
				if (iinod<jjnod)
				{
					inod = iinod;
					jnod = jjnod;
				}
				else
				{
					inod = jjnod;
					jnod = iinod;
				}
				xip = (xyz[0][inod-1] + xyz[0][jnod-1]) / 2.0;
				yip = (xyz[1][inod-1] + xyz[1][jnod-1]) / 2.0;
				zip = (xyz[2][inod-1] + xyz[2][jnod-1]) / 2.0;
				if (fabs(zip)<h)
				{
					for(na=0;na<nwperi;na++)
					{
						ja = iwperi[(int)na/100][(int)fmod(na,100)];
						ke = iwelem[(int)ja/100][(int)fmod(ja,100)];
						le = iwedge[(int)ja/100][(int)fmod(ja,100)];
						kknod = nod[inode1[le][0]][ke];
						llnod = nod[inode1[le][0]][ke];
						if (kknod<llnod)
						{
							knod = kknod;
							lnod = llnod;
						}
						else
						{
							knod = llnod;
							lnod = kknod;
						}
						xjp = (xyz[0][knod-1] + xyz[0][lnod-1]) / 2.0;
						yjp = (xyz[1][knod-1] + xyz[1][lnod-1]) / 2.0;
						zjp = (xyz[2][knod-1] + xyz[2][lnod-1]) / 2.0;
						dxy = sqrt((xip - xjp)*(xip-xjp) + (yip - yjp)*(yip-yjp));
						dz = fabs(zip - zjp);
						if ((dxy<h)&&(dz>=h))
						{
							iperi[0][np]=m;
							iperi[1][np]=ja;
							np=np+1;
						}
					}
				}
			}
		}
        if (jk==0)
		{
			nwperi=0;
			ne=*nedge;
			for(n=0;n<ne-1;n++)
			{
				i = iwelem[(int)n/100][(int)fmod(n,100)];
				j = iwedge[(int)n/100][(int)fmod(n,100)];
            
				iinod = nod[inode1[j][0]][i];
				jjnod = nod[inode1[j][1]][i];
				if (iinod<jjnod)
				{
					inod = iinod;
					jnod = jjnod;
				}
				else
				{
					inod = jjnod;
					jnod = iinod;

				}
				dz=fabs(xyz[2][inod-1]-xyz[2][jnod-1]);
				if (dz>=h)
				{
					inb[(int)nb/100][(int)fmod(nb,100)] = n;
					anb[(int)nb/100][(int)fmod(nb,100)]=0.0;
					nb = nb + 1;
				}
				else
				{
					iwperi[(int)nwperi/100][(int)fmod(nwperi,100)]=n;
					nwperi=nwperi+1;
				}
			}
			for(n=0;n<nwperi;n++)
			{
				m=iwperi[(int)n/100][(int)fmod(n,100)];
			}

			for(n=0;n<nwperi;n++)
			{
				m=iwperi[(int)n/100][(int)fmod(n,100)];
				i=iwelem[(int)m/100][(int)fmod(m,100)];
				j=iwedge[(int)m/100][(int)fmod(m,100)];
				iinod = nod[inode1[j][0]][i];
				jjnod = nod[inode1[j][1]][i];
				if (iinod<jjnod)
				{
					inod = iinod;
					jnod = jjnod;
				}
				else
				{
					inod = jjnod;
					jnod = iinod;
				}
				xip = (xyz[0][inod-1]+xyz[0][jnod-1]) / 2.0;
				yip = (xyz[1][inod-1]+xyz[1][jnod-1]) / 2.0;
				zip = (xyz[2][inod-1]+xyz[2][jnod-1]) / 2.0;
				if (fabs(zip)<h)
				{
					for(na=0;na<nwperi;na++)
					{
						ja = iwperi[(int)na/100][(int)fmod(na,100)];
						ke = iwelem[(int)ja/100][(int)fmod(ja,100)];
						le = iwedge[(int)ja/100][(int)fmod(ja,100)];
						kknod = nod[inode1[le][0]][ke];
						llnod = nod[inode1[le][0]][ke];
						if (kknod<llnod)
						{
							knod = kknod;
							lnod = llnod;
						}
						else
						{
							knod = llnod;
							lnod = kknod;
						}
						xjp = (xyz[0][knod-1] + xyz[0][lnod-1]) / 2.0;
						yjp = (xyz[1][knod-1] + xyz[1][lnod-1]) / 2.0;
						zjp = (xyz[2][knod-1] + xyz[2][lnod-1]) / 2.0;
						dxy = sqrt((xip - xjp)*(xip-xjp) + (yip - yjp)*(yip-yjp));
						dz = fabs(zip - zjp);
						if ((dxy<h)&&(dz>=h))
						{
							iperi[0][np]=m;
							iperi[1][np]=ja;
							np=np+1;
						}
					}
				}
			}
		}
		if(jk==8)
		{
			int jikux=ir[7];
			int jikuy=ir[8];
			int jikuz=ir[9];
			xmin=adiri[k][1]/fact-h;
			ymin=adiri[k][2]/fact-h;
			zmin=adiri[k][3]/fact-h;
			xmax=adiri[k][1]/fact+h;
			ymax=adiri[k][2]/fact+h;
			zmax=adiri[k][3]/fact+h;
			for(n=0;n<*nedge;n++)
			{
				int ie=iwelem[(int)n/100][(int)fmod(n,100)];
				int je=iwedge[(int)n/100][(int)fmod(n,100)];
				iinod=nod[inode1[je][0]][ie];
				jjnod=nod[inode1[je][1]][ie];
				if(iinod<jjnod)
				{
					inod=iinod;
					jnod=jjnod;
				}
				else
				{
					inod=jjnod;
					jnod=iinod;
				}
				xip=(xyz[0][inod-1] + xyz[0][jnod-1]) / 2.0;
				yip=(xyz[1][inod-1] + xyz[1][jnod-1]) / 2.0;
				zip=(xyz[2][inod-1] + xyz[2][jnod-1]) / 2.0;
				if(((xip>xmin)&&(xip<xmax))&&((yip>ymin)&&(yip<ymax))&&(zp>0))
				{
					inb[(int)nb/100][(int)fmod(nb,100)]=n;
					anb[(int)nb/100][(int)fmod(nb,100)]=fdiri(ie,je,adir,nod,xyz,npoint,nelem,nvert,ndim,nedges);
					nb=nb+1;
					if(iinod>=jjnod)
					{
						anb[(int)nb/100][(int)fmod(nb,100)]=-anb[(int)nb/100][(int)fmod(nb,100)];
					}
				}
			}
		}
		if(jk==9)
		{
			h=5.0;
			zahyo=(adiri[k][0]/fact)*(adiri[k][0]/fact);
			adir[0]=adiri[k][1];
			adir[1]=adiri[k][2];
			adir[2]=adiri[k][3];
			for(n=0;n<*nedge;n++)
			{
				int ie=iwelem[(int)n/100][(int)fmod(n,100)];
				int je=iwedge[(int)n/100][(int)fmod(n,100)];
				iinod=nod[inode1[je][0]][ie];
				jjnod=nod[inode1[je][1]][ie];
				if(iinod<jjnod)
				{
					inod=iinod;
					jnod=jjnod;
				}
				else
				{
					inod=jjnod;
					jnod=iinod;
				}
				rpi=xyz[0][inod-1]*xyz[0][inod-1]+xyz[1][inod-1]*xyz[1][inod-1];
				rpj=xyz[0][jnod-1]*xyz[0][jnod-1]+xyz[1][jnod-1]*xyz[1][jnod-1];
			}
		}
		idiri[k][1]=nb;
	}
	*nboun=nb;
	*nperi=np;
	free(adir);
	return 0;
}


int initao(int *ifsin,double **ao,int *nedge)
{
	int i;

	for(i=0;i<=*nedge;i++)
	{
		ao[(int)i/100][(int)fmod(i,100)]=0.0;
	}

	return 0;
}

double fdiri(int i,int j,double *adir,int **nod,double **xyz,int npoint,int nelem,int *nvert,int *ndim,int *nedges)
{
	double *u,*x,**dx,**dxi,**rn,**xn;
	double fj,fd,ad1,ad2,ad3;
	int jx,ii;

	u=(double *)malloc(sizeof(double)*3);
	x=(double *)malloc(sizeof(double)*3);
	dx=(double **)malloc(sizeof(double*)*3);
	for(ii=0;ii<3;ii++)
	{
		dx[ii]=(double *)malloc(sizeof(double)*3);
	}
	dxi=(double **)malloc(sizeof(double*)*3);
	for(ii=0;ii<3;ii++)
	{
		dxi[ii]=(double *)malloc(sizeof(double)*3);
	}
	xn=(double **)malloc(sizeof(double*)*36);
	for(ii=0;ii<36;ii++)
	{
		xn[ii]=(double *)malloc(sizeof(double)*3);
	}
	rn=(double **)malloc(sizeof(double*)*36);
	for(ii=0;ii<36;ii++)
	{
		rn[ii]=(double *)malloc(sizeof(double)*3);
	}

	u[0]=coefe1[j][0];
	u[1]=coefe1[j][1];
	u[2]=coefe1[j][2];

	fj=fjacob(0,2,i,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);

	free(u);
	free(x);
	for(ii=0;ii<3;ii++)
	{
		free(dxi[ii]);
	}
	free(dxi);
	for(ii=0;ii<36;ii++)
	{
		free(xn[ii]);
	}
	free(xn);
	for(ii=0;ii<36;ii++)
	{
		free(rn[ii]);
	}
	free(rn);

	jx=ixyze1[j];
	ad1=adir[0];
	ad2=adir[1];
	ad3=adir[2];

	fd=2.0*(ad1*dx[0][jx]+ad2*dx[1][jx]+ad3*dx[2][jx]);

	for(ii=0;ii<3;ii++)
	{
		free(dx[ii]);
	}
	free(dx);

	return fd;
}

double fjacob(int idebug,int itype,int i,double *u,double *x,double **dx,		\
			  double **dxi,double **xn,double **rn,int **nod,double **xyz,		\
			  int npoint,int nelem,int *nvert,int *ndim,int *nedges)
{
	int nd,nv,ndd,id,jd,je,jf,ned;
	int jx1,jx2,ix0,ix1,ix2;
	double sum,fj,fa,fb;

	if((int)fmod((int)itype/pow(2,0),2)==1)
	{
		for(nd=0;nd<*ndim;nd++)
		{
			sum=0.0;
			for(nv=0;nv<*nvert;nv++)
			{
				ndd=nod[nv][i];
				sum=sum+xyz[nd][ndd-1]*(1.0+coefn1[nv][0]*u[0])*(1.0+coefn1[nv][1]*u[1])*(1.0+coefn1[nv][2]*u[2]);
			}
			x[nd]=sum/8.0;
		}
	}
	if((int)fmod((int)itype/pow(2,1),2)==1)
	{
		for(id=0;id<*ndim;id++)
		{
			for(jd=0;jd<*ndim;jd++)
			{
				jx1=ir[jd+1];
				jx2=ir[jd+2];
				sum=0.0;
				for(nv=0;nv<*nvert;nv++)
				{
					ndd=nod[nv][i];
					sum=sum+(1.0+coefn1[nv][jx1]*u[jx1])*(1.0+coefn1[nv][jx2]*u[jx2])*coefn1[nv][jd]*xyz[id][ndd-1];
				}
				dx[id][jd]=sum/8.0;
			}
		}
		fj=dx[0][0]*(dx[1][1]*dx[2][2]-dx[1][2]*dx[2][1])+dx[0][1]*(dx[1][2]*dx[2][0]-dx[1][0]*dx[2][2])+dx[0][2]*(dx[1][0]*dx[2][1]-dx[1][1]*dx[2][0]);
		
		if(fj<=0.0)
		{
			printf("fjacob(%d)=%e\n",i,fj);
			printf("ERRER! fjacob lower than 0.0!\n");
			return 0;
		}
	}

	if((int)fmod((int)itype/pow(2,2),2)==1)
	{
		for(id=0;id<*ndim;id++)
		{
			ix1=ir[id+1];
			ix2=ir[id+2];
			for(jd=0;jd<*ndim;jd++)
			{
				jx1=ir[jd+1];
				jx2=ir[jd+2];
				dxi[id][jd]=(dx[jx1][ix1]*dx[jx2][ix2]-dx[jx1][ix2]*dx[jx2][ix1])/fj;
			}
		}
	}
	if((int)fmod((int)itype/pow(2,3),2)==1)
	{
		for(je=0;je<*nedges;je++)
		{
			ix0=ixyze1[je];
			ix1=ir[ix0+1];
			ix2=ir[ix0+2];
			for(nd=0;nd<*ndim;nd++)
			{
				xn[je][nd]=dxi[ix0][nd]*(1.0+coefe1[je][ix1]*u[ix1])*(1.0+coefe1[je][ix2]*u[jx2])/8.0;
			}
		}
	}
	if((int)fmod((int)itype/pow(2,4),2)==1)
	{
		for(jf=0;jf<*nedges;jf++)
		{
			ix0=ixyze1[jf];
			ix1=ir[ix0+1];
			ix2=ir[ix0+2];
			for(ned=0;ned<*ndim;ned++)
			{
				jx1=ir[ned+1];
				jx2=ir[ned+2];
				fa=dxi[ix0][jx2]*dxi[ix1][jx1]-dxi[ix0][jx1]*dxi[ix1][jx2];
				fb=dxi[ix0][jx2]*dxi[ix2][jx1]-dxi[ix0][jx1]*dxi[ix2][jx2];
				rn[jf][ned]=(coefe1[jf][ix1]*(1.0+coefe1[jf][ix2]*u[ix2])*fa+coefe1[jf][ix2]*(1.0+coefe1[jf][ix1]*u[ix1])*fb)/8.0;
			}
		}
	}
	
	return fj;
}

int fem(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,int *nedges,			\
		int **idiri,int *ndiri,int **iperi,int *nperi,int **inb,double **anb,int *nboun,int **imate,		\
		int *nblok,double **amate,int *nmate,int **icur,int *ncur,int *nmax,								\
		double *dmin,int *icga,int *nicmax,double *dicmin,double *freq,int *ifsin,int *nstep,				\
		double *dt,double **ao,int *nedge,double *xi,double *w,int ngauss,									\
		int ntmax,int **iwh,int **iwcic,int **iwg,int **iwgu,int nsum,int ndebug,int magmax,int *nmag,		\
		int **imag,double **amag,int *ndirif,int **idirif,int *nusef)
{
	int **ifa,**ifaf,**inbf;
	int ned,neds,nv,nd,**nk,**nkp;
	int i,istep,**nump,iihn,ibl,j;
	int iss,iee,il,ntp,nband;
	int nbounf,nficcg,ifcon,nnr;
	int nfeddy,nfpot,nfflux;
	double ammesh,ampot,amiccg,amtotal;
	double ts;
	double **fdb,afsin;
	double **tpot,**ttpot,**anbf,**a;
	double **h,**cic,**g,**gu;
	double **aa,*bo;
	double **p,**r,**s,**ai;
	double t0,t1,t2,t3,t4,t5,tg,tic,te;
	FILE *ficcg,*fcur,*fpot,*fflux,*feddy,*fa;

	printf("\n\n program FEM start...\n\n");
	t0=gettime();
	ned=*nedge;
	neds=*nedges;
	nv=*nvert;
	nd=*ndim;
	nficcg=1;
	nfeddy=2;
	nfpot=3;
	nfflux=4;
	if(nficcg>0)
	{
		ficcg=fopen("cpuiccg.dat","w");
	}

	ifa=(int**)malloc(sizeof(int*)*((int)npoint/25+1));
	for(i=0;i<(int)npoint/25+1;i++)
	{
		ifa[i]=(int*)malloc(sizeof(int)*100);
	}
	inbf=(int**)malloc(sizeof(int*)*((int)npoint/100));
	for(i=0;i<(int)npoint/100;i++)
	{
		inbf[i]=(int*)malloc(sizeof(int)*100);
	}
	ifaf=(int**)malloc(sizeof(int*)*((int)npoint/100));
	for(i=0;i<(int)npoint/100;i++)
	{
		ifaf[i]=(int*)malloc(sizeof(int)*100);
	}

	ptifa(idiri,ndiri,iperi,inb,nboun,ifa,nedge,ntmax,
		  &ntp,ndebug,npoint,ifaf,&nbounf,inbf,nblok,
		  nelem,nvert,imate,nod,nusef,nficcg,ficcg,nperi);

	printf("subroutine ptifa ok!\n");

	if(*icga==0)
	{
		bandcl();
	}
	else if(*icga==1)
	{
		nk=(int**)malloc(sizeof(int*)*((int)npoint/25+1));
		for(i=0;i<(int)npoint/25+1;i++)
		{
			nk[i]=(int*)malloc(sizeof(int)*100);
		}
		nkp=(int**)malloc(sizeof(int*)*((int)npoint/25+1));
		for(i=0;i<(int)npoint/25+1;i++)
		{
			nkp[i]=(int*)malloc(sizeof(int)*100);
		}
		nump=(int**)malloc(sizeof(int*)*((int)npoint/25*19));
		for(i=0;i<(int)npoint/25*19;i++)
		{
			nump[i]=(int*)malloc(sizeof(int)*100);
		}
		iccgnb(nelem,edge,nedges,ifa,nedge,nk,nkp,ntmax,
			   &ntp,nump,&nband,iwh,iwcic,iwg,iwgu,nsum,
			   ndebug,nficcg,nvert,nod,npoint,ifaf,nusef,ficcg);

		printf("subroutine iccgnb ok!\n");
		ammesh=(nelem*(nv+neds)*4+npoint*nd*8)/(1024*1024);
		
		if(*nusef==0)
		{
			ampot=((4+8*2)*ned)/(1024*1024);
		}
		else if(*nusef==1)
		{
			ampot=((4+8*2)*(ned+npoint))/(1024*1024);
		}
		amiccg=((4+8*2)*nband+(4*2+8*5)*ntp)/(1024*1024);
		amtotal=ammesh+ampot+amiccg;

		printf("mesh(nod,xyz,edge)\t\t%e\n",ammesh);
		printf("pot(ifa,a,ao)\t\t\t%e\n",ampot);
		printf("iccg(nump,h,cic,nk,nkp,g,aic,sic,ric,pic)\t%e\n",amiccg);
		printf("total(mesh,pot,iccg)\t\t%e\n",amtotal);

		fprintf(ficcg,"mesh(nod,xyz,edge)\t\t%e\n",ammesh);
		fprintf(ficcg,"pot(ifa,a,ao)\t\t\t%e\n",ampot);
		fprintf(ficcg,"iccg(nump,h,cic,nk,nkp,g,aic,sic,ric,pic)\t%e\n",amiccg);
		fprintf(ficcg,"total(mesh,pot,iccg)\t\t%e\n",amtotal);
	}
	if(*ncur>0)
	{
		tpot=(double **)malloc(sizeof(double*)*((int)(npoint/25+1)));
		for(i=0;i<(int)(npoint/25+1);i++)
		{
			tpot[i]=(double *)malloc(sizeof(double)*100);
		}
		fcur=fopen("tpot.dat","r");
		for(i=0;i<*nedge/3+1;i++)
		{
			fscanf(fcur,"%le\t%le\t%le",&tpot[(int)3*i/100][(int)fmod(3*i,100)],&tpot[(int)(3*i+1)/100][(int)fmod(3*i+1,100)],&tpot[(int)(3*i+2)/100][(int)fmod(3*i+2,100)]);
		}
		fclose(fcur);
	}
	fdb=(double**)malloc(sizeof(double*)*nelem);
	for(i=0;i<nelem;i++)
	{
		fdb[i]=(double*)malloc(sizeof(double)*3);
	}
	anbf=(double**)malloc(sizeof(double*)*((int)npoint/100));
	for(i=0;i<(int)npoint/100;i++)
	{
		anbf[i]=(double*)malloc(sizeof(double)*100);
	}
	ttpot=(double**)malloc(sizeof(double*)*((int)npoint/25+1));
	for(i=0;i<(int)npoint/25+1;i++)
	{
		ttpot[i]=(double*)malloc(sizeof(double)*100);
	}
	gu=(double**)malloc(sizeof(double*)*((int)npoint/25+1));
	for(i=0;i<(int)npoint/25+1;i++)
	{
		gu[i]=(double*)malloc(sizeof(double)*100);
	}
	for(i=0;i<(int)npoint/25+1;i++)
	{
		for(j=0;j<100;j++)
		{
			gu[i][j]=0.0;
		}
	}
	cic=(double**)malloc(sizeof(double*)*((int)npoint/25*19));
	for(i=0;i<(int)npoint/25*19;i++)
	{
		cic[i]=(double*)malloc(sizeof(double)*100);
	}
	h=(double**)malloc(sizeof(double*)*((int)npoint/25*19));
	for(i=0;i<(int)npoint/25*19;i++)
	{
		h[i]=(double*)malloc(sizeof(double)*100);
	}
	g=(double**)malloc(sizeof(double*)*((int)npoint/25+1));
	for(i=0;i<(int)npoint/25+1;i++)
	{
		g[i]=(double*)malloc(sizeof(double)*100);
	}
	for(i=0;i<(int)npoint/25+1;i++)
	{
		for(j=0;j<100;j++)
		{
			g[i][j]=0.0;
		}
	}
	a=(double**)malloc(sizeof(double*)*((int)npoint/25+1));
	for(i=0;i<(int)npoint/25+1;i++)
	{
		a[i]=(double*)malloc(sizeof(double)*100);
	}
	for(i=0;i<(int)npoint/25+1;i++)
	{
		for(j=0;j<100;j++)
		{
			a[i][j]=0.0;
		}
	}
	aa=(double**)malloc(sizeof(double*)*((int)npoint/25+1));
	for(i=0;i<(int)npoint/25+1;i++)
	{
		aa[i]=(double*)malloc(sizeof(double)*100);
	}
	for(i=0;i<(int)npoint/25+1;i++)
	{
		for(j=0;j<100;j++)
		{
			aa[i][j]=0.0;
		}
	}
	bo=(double*)malloc(sizeof(double)*nelem);
	p=(double **)malloc(sizeof(double*)*(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		p[i]=(double *)malloc(sizeof(double)*100);
	}
	r=(double **)malloc(sizeof(double*)*(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		r[i]=(double *)malloc(sizeof(double)*100);
	}
	s=(double **)malloc(sizeof(double*)*(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		s[i]=(double *)malloc(sizeof(double)*100);
	}
	ai=(double **)malloc(sizeof(double*)*(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		ai[i]=(double *)malloc(sizeof(double)*100);
	}
	
	if(nfflux>0)
	{
		fflux=fopen("flux.dat","w");
	}
	
	if(nfpot>0)
	{
		fpot=fopen("pot.dat","w");
	}
	if(nfeddy>0)
	{
		feddy=fopen("eddy.dat","w");
	}
	fa=fopen("iccgtime.dat","w");
	fcur=fopen("iccg.dat","w");
	t1=gettime();
	fprintf(fa,"%e\n",(t1-t0));
	/*�X�e�b�v���[�v�X�^�[�g*/
	for(istep=1;istep<=*nstep;istep++)
	{
		tg=0.0;
		tic=0.0;
		te=0.0;
		t2=gettime();
		initst(nelem,ndim,fdb);

		printf("subroutine initst[%d] ok!\n",istep-1);

		wtcal(freq,ifsin,dt,istep,&afsin);

		printf("subroutine wtcal[%d] ok!\n",istep-1);


		setdir(idiri,ndiri,inb,anb,nboun,a,nedge,
			   ndebug,&afsin,ndirif,idirif,&nbounf,
			   inbf,anbf,npoint,nusef);

		printf("subroutine setdir[%d] ok!\n",istep-1);
/*�d�����͋y�ѓd�����ɂ����}�g���N�X�`��*/
		if(*ncur>0)
		{
			setcur(tpot,ttpot,nedges,nedge,nelem,ncur,&afsin);

			printf("subroutine setcur[%d] ok!\n",istep-1);


			gujcal(nod,nelem,nvert,xyz,npoint,ndim,edge,
				   nedges,icur,ncur,ifa,nedge,gu,
				   xi,w,ngauss,iwh,cic,nsum,ndebug,ttpot,&ntp);

			printf("subroutine gujcal[%d] ok!\n",istep-1);

		}
/*���΍��ɂ����}�g���N�X�`��*/
		if(*nmag>0)
		{
			gumcal(nod,nelem,nvert,xyz,npoint,ndim,edge,
				   nedges,nmag,imag,amag,ifa,nedge,gu,
				   xi,w,ngauss,iwh,cic,nsum,ndebug,nk,nkp,
				   nump,&nband,h,&ntp);

			printf("subroutine gumcal[%d] ok!\n",istep-1);

		}
/*�j���[�g������t�\���@���[�v�X�^�[�g*/
		for(iihn=1;iihn<=*nmax;iihn++)
		{
			initnr(&ntp,&nband,h,g,gu);
			if(*icga==0)
			{
				gss();

				printf("subroutine gss[%d][%d] ok!\n",istep-1,iihn-1);

			}
			else if(*icga==1)
			{
/*�W���}�g���N�X�̍쐬*/
				for(ibl=1;ibl<=*nblok;ibl++)
				{
					iss=imate[ibl-1][0];
					iee=imate[ibl-1][1];
					il=imate[ibl-1][2];
					if(imate[ibl-1][3]>0)
					{
						icghfi(nod,nelem,nvert,xyz,npoint,ndim,edge,nedges,
							   imate,nblok,amate,nmate,freq,ifa,a,ao,nedge,
							   xi,w,ngauss,nk,nkp,&ntp,nump,&nband,h,g,
							   ndebug,dt,&iss,&iee,&il,fcur);

						printf("subroutine icghfi[%d][%d][%d] ok!\n",istep-1,iihn-1,ibl-1);
					}

					icgcl(nod,nelem,nvert,xyz,npoint,ndim,edge,nedges,
						  imate,nblok,amate,nmate,ifa,a,nedge,xi,w,
						  ngauss,nk,nkp,&ntp,nump,&nband,h,g,ndebug,
						  iihn,&iss,&iee,&il,aa,istep);

					printf("subroutine icgcl[%d][%d][%d] ok!\n",istep-1,iihn-1,ibl-1);
				}
				t3=gettime();
				tg=tg+(t3-t2);
				iccg(&nband,&ntp,h,g,cic,nk,nkp,nump,nicmax,dicmin,
					 iihn,istep,nficcg,ficcg,npoint,s,ai,r,p);
				printf("subroutine iccg[%d][%d] ok!\n",istep-1,iihn-1);
				t4=gettime();
				tic=tic+(t4-t3);
			}
			banda(nod,nelem,nvert,xyz,npoint,ndim,edge,
				  nedges,ifa,a,nedge,fdb,bo,&ntp,g,aa,
				  iihn,xi,ifaf,nusef);

			printf("subroutine banda[%d][%d] ok!\n",istep-1,iihn-1);
			if(*nmax<=1)
			{
				goto istepend;
			}
			if(iihn<=1)
			{
				goto iihnend;
			}
			conver(nelem,ndim,dmin,fdb,bo,nficcg,iihn,&ifcon,ficcg);

			printf("subroutine conver[%d][%d] ok!\n",istep-1,iihn-1);

			if((iihn>1)&&(ifcon<=0))
			{
				printf("******convergence ok ! at step %d ******\n",istep);
				printf("------ ( newton-raphson method ) ------\n");
/*				fprintf(ficcg,"******convergence ok ! at step %d ******\n",istep);
				fprintf(ficcg,"------ ( newton-raphson method ) ------\n");*/
				goto istepend;
			}
iihnend:;
		}
		printf("******no convergence at step %d ******\n",istep);
		printf("------ ( newton-raphson method ) ------\n");
/*		fprintf(ficcg,"******no convergence at step %d ******\n",istep);
		fprintf(ficcg,"------ ( newton-raphson method ) ------\n");*/
istepend:;
		nnr=iihn;

		if(nfeddy>0)
		{
			eddy(nod,nelem,nvert,xyz,npoint,ndim,edge,nedges,
				 imate,nblok,amate,nmate,freq,a,ao,nedge,nfeddy,
					dt,nusef,feddy);
			printf("subroutine eddy[%d][%d] ok!\n",istep-1,iihn-1);
		}
		

		atoao(a,ao,nedge);

		printf("subroutine atoao[%d][%d] ok!\n",istep-1,iihn-1);
		if(nfflux>0)
		{
			for(i=1;i<=*nedge;i++)
			{
				fprintf(fpot,"%e\n",a[(int)i/100][(int)fmod(i,100)]);
			}
		}
		if(nfpot>0)
		{
			for(i=0;i<nelem;i++)
			{
				fprintf(fflux,"%e\t%e\t%e\n",fdb[i][0],fdb[i][1],fdb[i][2]);
			}
		}
		t5=gettime();
		te=te+(t5-t4);
		fprintf(fa,"%e\t%e\t%e\n",tg,tic,te);
	}
	fclose(fa);
	if(nfpot>0)
	{
		fclose(fpot);
	}
	if(nfflux>0)
	{
		fclose(fflux);
	}
	if(nfeddy>0)
	{
		fclose(feddy);
	}
	printf("file closed...\n");
	printf("ndiri=%d\n",*ndiri);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(ifa[i]);
	}
	free(ifa);
	for(i=0;i<(int)(npoint/100);i++)
	{
		free(ifaf[i]);
	}
	free(ifaf);
	for(i=0;i<(int)(npoint/100);i++)
	{
		free(inbf[i]);
	}
	free(inbf);
	for(i=0;i<(int)(npoint/100);i++)
	{
		free(anbf[i]);
	}
	free(anbf);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(a[i]);
	}
	free(a);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(aa[i]);
	}
	free(aa);
	if(*ncur>0)
	{
		for(i=0;i<(int)(npoint/25+1);i++)
		{
			free(tpot[i]);
		}
		free(tpot);
		for(i=0;i<(int)(npoint/25+1);i++)
		{
			free(ttpot[i]);
		}
		free(ttpot);
	}
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(nk[i]);
	}
	free(nk);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(nkp[i]);
	}
	free(nkp);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(g[i]);
	}
	free(g);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(gu[i]);
	}
	free(gu);
	for(i=0;i<(int)(npoint/25*19);i++)
	{
		free(nump[i]);
	}
	free(nump);
	for(i=0;i<(int)(npoint/25*19);i++)
	{
		free(h[i]);
	}
	free(h);
	for(i=0;i<(int)(npoint/25*19);i++)
	{
		free(cic[i]);
	}
	free(cic);
	for(i=0;i<nelem;i++)
	{
		free(fdb[i]);
	}
	free(fdb);
	free(bo);
	for(i=0;i<27;i++)
	{
		free(nod[i]);
	}
	free(nod);
	for(i=0;i<3;i++)
	{
		free(xyz[i]);
	}
	free(xyz);
	for(i=0;i<nelem;i++)
	{
		free(edge[i]);
	}
	free(edge);
	for(i=0;i<25;i++)
	{
		free(idiri[i]);
	}
	free(idiri);
	for(i=0;i<2;i++)
	{
		free(iperi[i]);
	}
	free(iperi);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(inb[i]);
	}
	free(inb);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(anb[i]);
	}
	free(anb);
	for(i=0;i<25;i++)
	{
		free(imate[i]);
	}
	free(imate);
	for(i=0;i<25;i++)
	{
		free(amate[i]);
	}
	free(amate);
	for(i=0;i<25;i++)
	{
		free(icur[i]);
	}
	free(icur);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(ao[i]);
	}
	free(ao);
	free(xi);
	free(w);
	fclose(ficcg);
	for(i=0;i<npoint/25+1;i++)
	{
		free(r[i]);
	}
	free(r);
	for(i=0;i<npoint/25+1;i++)
	{
		free(p[i]);
	}
	free(p);
	for(i=0;i<npoint/25+1;i++)
	{
		free(s[i]);
	}
	free(s);
	for(i=0;i<npoint/25+1;i++)
	{
		free(ai[i]);
	}
	free(ai);
	return 0;
}

int ptifa(int **idiri,int *ndiri,int **iperi,int **inb,int *nboun,int **ifa,int *nedge,int ntmax,
		  int *ntp,int ndebug,int npoint,int **ifaf,int *nbounf,int **inbf,int *nblok,
		  int nelem,int *nvert,int **imate,int **nod,int *nusef,int nficcg,FILE *ficcg,int *nperi)
{
	int ned,id,nper,i,nb,ifn;
	int ii,is,ie,j,k,n,ip1,ip2,ipe,ipeab;
	FILE *fa;

	for(ned=1;ned<=*nedge;ned++)
	{
		ifa[(int)ned/100][(int)fmod(ned,100)]=1;
	}

	for(id=0;id<*ndiri;id++)
	{
		for(nb=idiri[id][0];nb<idiri[id][1];nb++)
		{
			ifn=inb[(int)nb/100][(int)fmod(nb,100)]+1;
			ifa[(int)ifn/100][(int)fmod(ifn,100)]=0;
		}
	}
	for(nper=0;nper<*nperi;nper++)
	{
		ipe=iperi[1][nper]+1;
		ipeab=(int)fabs(ipe);
		ifa[(int)ipeab/100][(int)fmod(ipeab,100)]=0;
	}

	if(*nusef==1)
	{
		for(i=1;i<=npoint;i++)
		{
			ifaf[(int)i/100][(int)fmod(i,100)]=0;
		}
		for(i=0;i<*nblok;i++)
		{
			ii=imate[i][3];
			if(ii>0)
			{
				is=imate[i][0];
				ie=imate[i][1];
				for(j=is;j<=ie;j++)
				{
					for(k=0;k<*nvert;k++)
					{
						n=nod[k][j-1];
						ifaf[(int)n/100][(int)fmod(n,100)]=1;
					}
				}
			}
		}
	}
	*ntp=0;
	for(ned=1;ned<=*nedge;ned++)
	{
		if(ifa[(int)ned/100][(int)fmod(ned,100)]>0)
		{
			*ntp=*ntp+1;
			ifa[(int)ned/100][(int)fmod(ned,100)]=*ntp;
		}
	}
	if(*nusef==1)
	{
		for(ned=1;ned<=npoint;ned++)
		{
			if(ifaf[(int)ned/100][(int)fmod(ned,100)]>0)
			{
				*ntp=*ntp+1;
				ifaf[(int)ned/100][(int)fmod(ned,100)]=*ntp;
			}
		}
	}

	printf("ntp(matrix leng) = %d\n",*ntp);
	fprintf(ficcg,"ntp(matrix leng) = %d\n",*ntp);

	for(nper=0;nper<*nperi;nper++)
	{
		ip1=iperi[0][nper]+1;
		ip2=iperi[1][nper]+1;
		if(ip2>0)
		{
			ifa[(int)ip2/100][(int)fmod(ip2,100)]=ifa[(int)ip1/100][(int)fmod(ip1,100)];
		}
		if(ip2<0)
		{
			ifa[(int)ip2/100][(int)fmod(ip2,100)]=-ifa[(int)ip1/100][(int)fmod(ip1,100)];
		}
	}
	fa=fopen("ifa.dat","w");
	for(ned=1;ned<=*nedge;ned++)
	{
		fprintf(fa,"%d\t%d\n",ned,ifa[(int)ned/100][(int)fmod(ned,100)]);
	}
	fclose(fa);


	return 0;
}

int bandcl(void)
{
	return 0;
}

int iccgnb(int nelem,int **edge,int *nedges,int **ifa,int *nedge,int **nk,int **nkp,int ntmax,
		   int *ntp,int **nump,int *nband,int **iwork,int **iwelem,int **iwk,int **iwkp,int nsum,
		   int ndebug,int nficcg,int *nvert,int **nod,int npoint,int **ifaf,int *nusef,FILE *ficcg)
{
	int **iwflag;
	int i,nt,isum,ia,ja,ne,ip,is,ed,nd,ned;
	int nend,nkmax,istart,nstart,nb,nba;
	int ifaj,ifai,isort;
	int iskp,iswk,nup;
	FILE *fb;

	iwflag=(int**)malloc(sizeof(int*)*ntmax/100);
	for(i=0;i<ntmax/100;i++)
	{
		iwflag[i]=(int*)malloc(sizeof(int)*100);
	}

	for(nt=0;nt<=*ntp;nt++)
	{
		iwk[(int)nt/100][(int)fmod(nt,100)]=0;
	}
	for(i=1;i<=*nedge*2;i++)
	{
		iwork[(int)i/100][(int)fmod(i,100)]=0;
	}

	isum=0;
	ned=*nedges;
	for(ia=0;ia<nelem;ia++)
	{
		for(ja=0;ja<ned;ja++)
		{
			ed=edge[ia][ja]+1;
			ifai=ifa[(int)ed/100][(int)fmod(ed,100)];
			if((ifai>0)/*&&(ifai<*nedge)*/)
			{
				iwork[(int)(isum+1)/100][(int)fmod(isum+1,100)]=ifai;
				iwork[(int)(isum+2)/100][(int)fmod(isum+2,100)]=ia;
				iwk[(int)ifai/100][(int)fmod(ifai,100)]=iwk[(int)ifai/100][(int)fmod(ifai,100)]+1;
				isum=isum+2;
			}
		}
		if(*nusef==1)
		{
			for(ja=0;ja<*nvert;ja++)
			{
				nd=nod[ja][ia];
				ifai=ifaf[(int)nd/100][(int)fmod(nd,100)];
				if((ifai>0)/*&&(ifai<*nedge)*/)
				{
					iwork[(int)isum/100][(int)fmod(isum,100)]=ifai;
					iwork[(int)(isum+1)/100][(int)fmod(isum+1,100)]=ia;
					iwk[(int)ifai/100][(int)fmod(ifai,100)]=iwk[(int)ifai/100][(int)fmod(ifai,100)]+1;
					isum=isum+2;
				}
			}
		}
	}

	iwkp[0][1]=1;
	for(ip=2;ip<=*ntp;ip++)
	{
		iwkp[(int)ip/100][(int)fmod(ip,100)]=iwkp[(int)(ip-1)/100][(int)fmod(ip-1,100)]+iwk[(int)(ip-1)/100][(int)fmod(ip-1,100)];
	}

	for(is=1;is<=isum;)
	{
		iswk=iwork[(int)is/100][(int)fmod(is,100)];
		iskp=iwkp[(int)iswk/100][(int)fmod(iswk,100)];
		iwelem[(int)iskp/100][(int)fmod(iskp,100)]=iwork[(int)(is+1)/100][(int)fmod(is+1,100)]+1;
		iwkp[(int)iswk/100][(int)fmod(iswk,100)]=iwkp[(int)iswk/100][(int)fmod(iswk,100)]+1;
		is=is+2;
	}
	for(nt=1;nt<=*ntp;nt++)
	{
		iwflag[(int)nt/100][(int)fmod(nt,100)]=0;
	}
	fb=fopen("femiwk.dat","w");
	for(i=1;i<=isum;i++)
	{
		fprintf(fb,"%d\t%d\t%d\t",i,iwork[(int)i/100][(int)fmod(i,100)],iwelem[(int)i/100][(int)fmod(i,100)]);
		if(i<*ntp)
		{
			fprintf(fb,"%d\t%d\n",iwk[(int)i/100][(int)fmod(i,100)],iwkp[(int)i/100][(int)fmod(i,100)]);
		}
		else
		{
			fprintf(fb,"\n");
		}
	}
	fclose(fb);


	*nband=0;
	nb=*nband;
	nend=0;
	nkmax=0;
	for(ifai=1;ifai<=*ntp;ifai++)
	{
		nk[(int)ifai/100][(int)fmod(ifai,100)]=0;
		istart=*nband+1;
		nstart=nend+1;
		nend=nstart+iwk[(int)ifai/100][(int)fmod(ifai,100)]-1;
		for(ne=nstart;ne<=nend;ne++)
		{
			ia=iwelem[(int)ne/100][(int)fmod(ne,100)];
			for(ja=0;ja<*nedges;ja++)
			{
				ed=edge[ia-1][ja]+1;
				ifaj=ifa[(int)ed/100][(int)fmod(ed,100)];
				if((ifaj>0)&&(ifai>=ifaj)&&(iwflag[(int)ifaj/100][(int)fmod(ifaj,100)]==0))
				{
					nb=nb+1;
					nump[(int)(nb)/100][(int)fmod(nb,100)]=ifaj;
					nk[(int)ifai/100][(int)fmod(ifai,100)]=nk[(int)ifai/100][(int)fmod(ifai,100)]+1;
					iwflag[(int)ifaj/100][(int)fmod(ifaj,100)]=1;
				}
			}
/*			if(*nusef==1)
			{
				for(ja=0;ja<*nvert;ja++)
				{
					nd=nod[ja][ia];
					ifaj=ifaf[(int)(nd-1)/100][(int)fmod(nd-1,100)];
					if((ifaj>0)&&(ifai>=ifaj)&&(iwflag[(int)ifaj/100][(int)fmod(ifaj,100)]==0))
					{
						nb=nb+1;
						nump[(int)(nb-1)/100][(int)fmod(nb-1,100)]=ifaj;
						nk[(int)ifai/100][(int)fmod(ifai,100)]=nk[(int)ifai/100][(int)fmod(ifai,100)]+1;
						iwflag[(int)ifaj/100][(int)fmod(ifaj,100)]=1;
					}
				}
			}	*/
		}
		*nband=nb;
		isort=2;
		if(isort==1)
		{
			shell(nump,ntmax*19,istart,nb);
		}
		if(isort==2)
		{
			bubble(nump,ntmax*19,istart,nb);
		}

		if(nkmax<nk[(int)ifai/100][(int)fmod(ifai,100)])
		{
			nkmax=nk[(int)ifai/100][(int)fmod(ifai,100)];
		}

		for(nba=istart;nba<=*nband;nba++)
		{
			nup=nump[(int)nba/100][(int)fmod(nba,100)];
			iwflag[(int)nup/100][(int)fmod(nup,100)]=0;
		}
	}
	for(i=0;i<ntmax/100;i++)
	{
		free(iwflag[i]);
	}
	free(iwflag);
	nkp[0][0]=0;
	nkp[0][1]=1;
	for(i=2;i<=*nedge;i++)
	{
		nkp[(int)i/100][(int)fmod(i,100)]=nkp[(int)(i-1)/100][(int)fmod(i-1,100)]+nk[(int)(i-1)/100][(int)fmod(i-1,100)];
	}

	printf("*************** iccg method ***************\n");
	printf("nband(all band width) = %d\n",*nband);
	printf("nkmax(max band width) = %d\n",nkmax);

	fprintf(ficcg,"*************** iccg method ***************\n");
	fprintf(ficcg,"nband(all band width) = %d\n",*nband);
	fprintf(ficcg,"nkmax(max band width) = %d\n",nkmax);

	return 0;
}
/* �������x�Ɠd���̏����� */
int initst(int nelem,int *ndim,double **fdb)
{
	int ne,nd;

	for(ne=0;ne<nelem;ne++)
	{
		for(nd=0;nd<*ndim;nd++)
		{
			fdb[ne][nd]=0.0e0;
		}
	}

	return 0;
}
/* �d���g�`�̌��� */
int wtcal(double *freq,int *ifsin,double *dt,int istep,double *afsin)
{
	double pi,timewa,d,f;

	d=*dt;
	f=*freq;
	pi=acos(-1.0);
	timewa=d*istep;
/* ���� */
	if(*ifsin==0)
	{
		*afsin=1.0;
	}
/* �ߓn����(����) */
	if(*ifsin==1)
	{
		*afsin=exp(-timewa/f);
	}
/* �ߓn����(����) */
	if(*ifsin==2)
	{
		*afsin=1.0-exp(-timewa/f);
	}
/* �����g�` */
	if(*ifsin==3)
	{
		*afsin=sqrt(2.0)*sin(2.0*pi*f*timewa);
	}
/* �ꎟ�㏸ */	
	if(*ifsin==4)
	{
		*afsin=timewa;
	}

	printf("Now be in subroutine wtcal...\n");
	printf("ifsin = %d\n",*ifsin);
	printf("istep = %d\n",istep);
	printf("afsin = %e\n",*afsin);
	printf("freq  = %e\n",*freq);
	printf("timewa= %e\n",timewa);

	return 0;
}
/* ���E��ɂ����x�N�g���|�e���V�����̌��� */
int setdir(int **idiri,int *ndiri,int **inb,double **anb,int *nboun,double **a,int *nedge,
		   int ndebug,double *afsin,int *ndirif,int **idirif,int *nbounf,
		
		   int **inbf,double **anbf,int npoint,int *nusef)
{
	int ned,id,nb,ib;
	double afs;
	FILE *fa;

	afs=*afsin;

	for(ned=0;ned<=*nedge;ned++)
	{
		a[(int)ned/100][(int)fmod(ned,100)]=0.0;
	}

	for(id=0;id<*ndiri;id++)
	{
		for(nb=idiri[id][0];nb<idiri[id][1];nb++)
		{
			ib=inb[(int)nb/100][(int)fmod(nb,100)]+1;
			a[(int)ib/100][(int)fmod(ib,100)]=1.0*anb[(int)nb/100][(int)fmod(nb,100)];
		}
	}
	fa=fopen("setdir.dat","w");
		for(ned=0;ned<=*nedge;ned++)
	{
		fprintf(fa,"%d\t%e\n",ned,a[(int)ned/100][(int)fmod(ned,100)]);
	}
		fclose(fa);

	return 0;
}
/* �X�e�b�v���̓d���x�N�g���|�e���V�����̌v�Z */
int setcur(double **tpot,double **ttpot,int *nedges,int *nedge,int nelem,int *ncur,double *afsin)
{
	int je;
	double afs;

	afs=*afsin;

	for(je=0;je<*nedge;je++)
	{
		ttpot[(int)je/100][(int)fmod(je,100)]=1.0*tpot[(int)je/100][(int)fmod(je,100)];
	}

	return 0;
}
/* �d�����ɂ����W���}�g���N�X�ƉE�Ӄx�N�g���̍쐬 */
int gujcal(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,
		   int *nedges,int **icur,int *ncur,int **ifa,int *nedge,double **gu,
		   double *xi,double *w,int ngauss,int **lin,double **dataj0,int nsum,int ndebug,double **ttpot,int *ntp)
{
	int ic,ie,i,j,k,im,je,ke,ned2;
	int im2,nod1,nod2,ed,nt;
	int iedges;
	double ww,hh;
	double *u,*x,**dx,**dxi,**xn,**rn;
	iedges=*nedges;

	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}

	for(nt=1;nt<=*ntp;nt++)
	{
		gu[(int)nt/100][(int)fmod(nt,100)]=0.0;
	}
	u=(double*)malloc(sizeof(double)*3);

	for(ic=0;ic<*ncur;ic++)
	{
		for(ie=icur[ic][0];ie<=icur[ic][1];ie++)
		{
			for(i=0;i<ngauss;i++)
			{
				u[0]=xi[i];
				for(j=0;j<ngauss;j++)
				{
					u[1]=xi[j];
					for(k=0;k<ngauss;k++)
					{
						u[2]=xi[k];
						ww=w[i]*w[j]*w[k]*fjacob(0,31,ie-1,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
						for(je=0;je<*nedges;je++)
						{
							ed=edge[ie-1][je]+1;
							im=ifa[(int)ed/100][(int)fmod(ed,100)];
							if(im>0)
							{
								for(ke=0;ke<*nedges;ke++)
								{
									ned2=edge[ie-1][ke]+1;
									hh=ww*(xn[je][0]*rn[ke][0]+xn[je][1]*rn[ke][1]+xn[je][2]*rn[ke][2])*ttpot[(int)(ned2-1)/100][(int)fmod(ned2-1,100)];
									im2=ifa[(int)ned2/100][(int)fmod(ned2,100)];

									nod1=nod[inode11[je][0]][ie-1];
									nod2=nod[inode11[je][1]][ie-1];
									if(nod1>nod2)
									{
										hh=-hh;
									}
									nod1=nod[inode11[ke][0]][ie-1];
									nod2=nod[inode11[ke][1]][ie-1];
									if(nod1>nod2)
									{
										hh=-hh;
									}
									gu[(int)im/100][(int)fmod(im,100)]=gu[(int)im/100][(int)fmod(im,100)]+hh;
								}
							}
						}
					}
				}
			}
		}
	}
	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	free(rn);
	free(u);
	free(x);

	return 0;
}
/* ���΍��Ɋւ����W���}�g���N�X�ƉE�Ӄx�N�g���̍쐬 */
int gumcal(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,
		   int *nedges,int *nmag,int **imag,double **amag,int **ifa,int *nedge,double **gu,
		   double *xi,double *w,int ngauss,int **iwh,double **cic,int nsum,int ndebug,int **nk,int **nkp,
		   int **nump,int *nband,double **h,int *ntp)
{
	int nt,iss,iee,ie,i,j,k,je,im,ke;
	int nod1,nod2,ibl,ed,iedges;
	double pi,xmu0,ww,hh;
	double **xnu,**xn,*u,*x,**dx,**dxi;
	double **rn,*xm0;
	FILE *fa;
	iedges=*nedges;

	for(nt=1;nt<=*ntp;nt++)
	{
		gu[(int)nt/100][(int)fmod(nt,100)]=0.0;
	}

	pi=acos(-1.0);
	xmu0=4.0*pi*0.0000001;
	
	printf("m0 = %e\n",xmu0);
	
	xnu=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		xnu[i]=(double*)malloc(sizeof(double)*3);
	}
	xnu[0][0]=1.0/xmu0;
	xnu[1][1]=1.0/xmu0;
	xnu[2][2]=1.0/xmu0;
	xm0=(double*)malloc(sizeof(double)*3);
	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}

	for(ibl=0;ibl<*nmag;ibl++)
	{
		iss=imag[ibl][0];
		iee=imag[ibl][1];

		printf("iss = %d\n",iss);
		printf("iee = %d\n",iee);

		for(ie=iss;ie<=iee;ie++)
		{
			for(i=0;i<ngauss;i++)
			{
				u[0]=xi[i];
				for(j=0;j<ngauss;j++)
				{
					u[1]=xi[j];
					for(k=0;k<ngauss;k++)
					{
						u[2]=xi[k];

						ww=w[i]*w[j]*w[k]*fjacob(0,31,ie-1,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
						for(je=0;je<*nedges;je++)
						{
							ed=edge[ie-1][je];
							im=ifa[(int)(ed+1)/100][(int)fmod(ed+1,100)];
							if(im>0)
							{
								for(ke=0;ke<*nedges;ke++)
								{
									m0cal(nod,nelem,nvert,xyz,npoint,ndim,
										  amag,nmag,ibl,xm0,ie-1,je,ke);
									/*printf("xm0=(%e,%e,%e)\n",xm0[0],xm0[1],xm0[2]);*/
									hh=ww*(xnu[0][0]*xn[je][0]*rn[ke][0]*xm0[0]+xnu[1][1]*xn[je][1]*rn[ke][1]*xm0[1]+xnu[2][2]*xn[je][2]*rn[ke][2]*xm0[2]);
									/*printf("hh=%e\n",hh);*/
									nod1=nod[inode11[je][0]][ie-1];
									nod2=nod[inode11[je][1]][ie-1];
									if(nod1>nod2)
									{
										hh=-hh;
									}
									nod1=nod[inode11[ke][0]][ie-1];
									nod2=nod[inode11[ke][1]][ie-1];
									if(nod1>nod2)
									{
										hh=-hh;
									}
									gu[(int)im/100][(int)fmod(im,100)]+=hh;
								}
							}
						}
					}
				}
			}
		}
	}
	fa=fopen("gumcal.dat","w");
	for(i=1;i<=*ntp;i++)
	{
		fprintf(fa,"%d\t%e\t%e\n",i,gu[(int)i/100][(int)fmod(i,100)],h[(int)i/100][(int)fmod(i,100)]);
	}
	fclose(fa);
	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	free(rn);
	free(u);
	free(x);
	free(xm0);

	return 0;
}

int m0cal(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,
		  double **amag,int *nmag,int ibl,double *xm0,int ie,int je,int ke)
{
	int nnod1,nnod2,nod1,nod2;
	double dx,dy,dz;
	nnod1=nod[inode11[ke][0]][ie];
	nnod2=nod[inode11[ke][1]][ie];

	if(nnod1<nnod2)
	{
		nod1=nnod1;
		nod2=nnod2;
	}
	else
	{
		nod1=nnod2;
		nod2=nnod1;
	}
	dx=xyz[0][nod2-1]-xyz[0][nod1-1];
	dy=xyz[1][nod2-1]-xyz[1][nod1-1];
	dz=xyz[2][nod2-1]-xyz[2][nod1-1];
	xm0[0]=amag[ibl][0]*dx;
	xm0[1]=amag[ibl][1]*dy;
	xm0[2]=amag[ibl][2]*dz;

	return 0;
}

int initnr(int *ntp,int *nband,double **h,double **g,double **gu)
{
	int nb,nt;
	FILE *fa;

	for(nb=1;nb<=*nband;nb++)
	{
		h[(int)nb/100][(int)fmod(nb,100)]=0.0;
	}
	fa=fopen("initnr.dat","w");
	for(nt=1;nt<=*ntp;nt++)
	{
		g[(int)nt/100][(int)fmod(nt,100)]=gu[(int)nt/100][(int)fmod(nt,100)];
		fprintf(fa,"%d\t%e\n",nt,g[(int)nt/100][(int)fmod(nt,100)]);
	}
	fclose(fa);

	return 0;
}

int gss(void)
{
	return 0;
}

int icghfi(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,int *nedges,
		   int **imate,int *nblok,double **amate,int *nmate,double *freq,int **ifa,double **a,double **ao,int *nedge,
		   double *xi,double *w,int ngauss,int **nk,int **nkp,int *ntp,int **nump,int *nband,double **h,double **g,
		   int ndebug,double *dt,int *iss,int *iee,int *il,FILE *fcur)
{
	int ie,i,j,k,je,ke,ned,im1;
	int noko,nod1,nod2,ned2,im2,iedges;
	double pi,delta,ww,hh,ddt,fj;
	double *u,*x,**dx,**dxi,**xn,**rn;
	FILE *fa;
	iedges=*nedges;

	pi=acos(-1.0);
	ddt=*dt;
	delta=amate[*il-1][3]/ddt;
	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}

	for(ie=*iss;ie<=*iee;ie++)
	{
		for(i=0;i<ngauss;i++)
		{
			u[0]=xi[i];
			for(j=0;j<ngauss;j++)
			{
				u[1]=xi[j];
				for(k=0;k<ngauss;k++)
				{
					u[2]=xi[k];
					fj=fjacob(0,14,ie-1,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
					ww=w[i]*w[j]*w[k]*fj;
					fprintf(fcur,"%e\n",fj);
					for(je=0;je<*nedges;je++)
					{
						ned=edge[ie-1][je]+1;
						im1=ifa[(int)ned/100][(int)fmod(ned,100)];
						if(im1>0)
						{
							for(ke=0;ke<*nedges;ke++)
							{
								ned2=edge[ie-1][ke]+1;
								im2=ifa[(int)ned2/100][(int)fmod(ned2,100)];
								hh=ww*delta*(xn[je][0]*xn[ke][0]+xn[je][1]*xn[ke][1]+xn[je][2]*xn[ke][2]);
								nod1=nod[inode11[je][0]][ie-1];
								nod2=nod[inode11[je][1]][ie-1];
								if(nod1>nod2)
								{
									hh=-hh;
								}
								nod1=nod[inode11[ke][0]][ie-1];
								nod2=nod[inode11[ke][1]][ie-1];
								if(nod1>nod2)
								{
									hh=-hh;
								}
								if((im2>0)&&(im2<=im1))
								{
									for(noko=nkp[(int)im1/100][(int)fmod(im1,100)];noko<nkp[(int)im1/100][(int)fmod(im1,100)]+nk[(int)im1/100][(int)fmod(im1,100)];noko++)
									{
										if(nump[(int)noko/100][(int)fmod(noko,100)]==im2)
										{
											h[(int)noko/100][(int)fmod(noko,100)]=h[(int)noko/100][(int)fmod(noko,100)]+hh;
											goto nodconti;
										}
									}
								}
nodconti:;
								g[(int)im1/100][(int)fmod(im1,100)]=g[(int)im1/100][(int)fmod(im1,100)]-hh*(a[(int)ned2/100][(int)fmod(ned2,100)]-ao[(int)ned2/100][(int)fmod(ned2,100)]);
							}
						}
					}
				}
			}
		}
	}
	fa=fopen("icghfi.dat","w");
	for(i=1;i<=*ntp;i++)
	{
		fprintf(fa,"%d\t%e\t%e\n",i,g[(int)i/100][(int)fmod(i,100)],h[(int)i/100][(int)fmod(i,100)]);
	}
	fclose(fa);
	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	free(rn);
	free(u);
	free(x);

	return 0;
}

int icgcl(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,int *nedges,
		  int **imate,int *nblok,double **amate,int *nmate,int **ifa,double **a,int *nedge,double *xi,double *w,
		  int ngauss,int **nk,int **nkp,int *ntp,int **nump,int *nband,double **h,double **g,int ndebug,
		  int iihn,int *iss,int *iee,int *il,double **aa,int istep)
{
	int melem,mgauss,ie,i,j,k,je,ke,le;
	int ned,ned2,im1,im2,ned3,nod1,nod2;
	int noko,il1,iedges;
	double pi,xmu0,ww,hh,tbb,cc,tcc;
	double ccc,ddx,xnux,bb;
	double *anyu,*dife;
	double *u,*x,**dx,**dxi,**xn,**rn,**xnu;
	FILE *fa,*fb;

	iedges=*nedges;

	pi=acos(-1.0);
	xmu0=4.0*pi*1.0e-7;

	melem=nelem;
	mgauss=ngauss;

	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	anyu=(double*)malloc(sizeof(double)*ngauss*ngauss*ngauss);
	dife=(double*)malloc(sizeof(double)*ngauss*ngauss*ngauss);
	for(i=0;i<3;i++)
	{
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}
	xnu=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		xnu[i]=(double*)malloc(sizeof(double)*3);
	}
	fb=fopen("rotn.dat","w");
	for(ie=*iss;ie<=*iee;ie++)
	{
		if((iihn>1)&&(*il==4))
		{
			nonline(melem,ie,ndim,a,edge,nedge,nedges,anyu,dife,
				    ngauss,xi,nod,xyz,npoint,nvert,aa);
		}
		else
		{
			il1=*il-1;
			xnu[0][0]=1.0/(xmu0*amate[il1][0]);
			xnu[1][1]=1.0/(xmu0*amate[il1][1]);
			xnu[2][2]=1.0/(xmu0*amate[il1][2]);
			fprintf(fb,"%e\t%e\t%e\n",xnu[0][0],xnu[1][1],xnu[2][2]);
		}
		for(i=0;i<ngauss;i++)
		{
			u[0]=xi[i];
			for(j=0;j<ngauss;j++)
			{
				u[1]=xi[j];
				for(k=0;k<ngauss;k++)
				{
					u[2]=xi[k];
					xnux=anyu[4*i+2*j+k];
					ddx=dife[4*i+2*j+k];
					ww=w[i]*w[j]*w[k]*fjacob(0,22,ie-1,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
					for(je=0;je<*nedges;je++)
					{
						ned=edge[ie-1][je]+1;
						im1=ifa[(int)ned/100][(int)fmod(ned,100)];
						if(im1>0)
						{
							for(ke=0;ke<*nedges;ke++)
							{
								ned2=edge[ie-1][ke]+1;
								im2=ifa[(int)ned2/100][(int)fmod(ned2,100)];
								if((iihn<=1)||(*il!=4))
								{
									hh=ww*(xnu[0][0]*rn[je][0]*rn[ke][0]+xnu[1][1]*rn[je][1]*rn[ke][1]+xnu[2][2]*rn[je][2]*rn[ke][2]);
									tbb=0.0;
								}
								else
								{
									tcc=0.0;
									tbb=0.0;
									hh=ww*(xnux*rn[je][0]*rn[ke][0]+xnux*rn[je][1]*rn[ke][1]+xnux*rn[je][2]*rn[ke][2]);
									for(le=0;le<*nedges;le++)
									{
										ned3=edge[ie-1][le]+1;
										cc=(rn[le][0]*rn[ke][0]+rn[le][1]*rn[ke][1]+rn[le][2]*rn[ke][2])*a[(int)ned3/100][(int)fmod(ned3,100)];
										nod1=nod[inode11[le][0]][ie-1];
										nod2=nod[inode11[le][1]][ie-1];
										if(nod1>nod2)
										{
											cc=-cc;
										}
										tcc=tcc+cc;
									}
									ccc=2.0*ddx*tcc;

									for(le=0;le<*nedges;le++)
									{
										ned3=edge[ie-1][le]+1;
										nod1=nod[inode11[le][0]][ie-1];
										nod2=nod[inode11[le][1]][ie-1];
										bb=ww*(rn[je][0]*rn[le][0]+rn[je][1]*rn[le][1]+rn[je][2]*rn[le][2])*ccc*a[(int)ned3/100][(int)fmod(ned3,100)];
										if(nod1>nod2)
										{
											bb=-bb;
										}
										tbb=tbb+bb;
									}
								}
								nod1=nod[inode11[je][0]][ie-1];
								nod2=nod[inode11[je][1]][ie-1];
								if(nod1>nod2)
								{
									hh=-hh;
									tbb=-tbb;
								}
								nod1=nod[inode11[ke][0]][ie-1];
								nod2=nod[inode11[ke][1]][ie-1];
								if(nod1>nod2)
								{
									hh=-hh;
									tbb=-tbb;
								}
								if((im2>0)&&(im2<=im1))
								{
									for(noko=nkp[(int)im1/100][(int)fmod(im1,100)];noko<nkp[(int)im1/100][(int)fmod(im1,100)]+nk[(int)im1/100][(int)fmod(im1,100)];noko++)
									{
										if(nump[(int)noko/100][(int)fmod(noko,100)]==im2)
										{
											h[(int)noko/100][(int)fmod(noko,100)]=h[(int)noko/100][(int)fmod(noko,100)]+hh+tbb;
											goto icgconti;
										}
									}
								}
icgconti:;
								fprintf(fb,"%d\t%d\t%d\t%d\t%d\t%e\t%e\n",i,j,k,je,ke,ww,hh);
								g[(int)im1/100][(int)fmod(im1,100)]=g[(int)im1/100][(int)fmod(im1,100)]-hh*a[(int)ned2/100][(int)fmod(ned2,100)];
							}
						}
					}
				}
			}
		}
	}
	fclose(fb);
	if(istep==1)
	{
		fa=fopen("iccghg.dat","w");
		for(i=1;i<=*nband;i++)
		{
			fprintf(fa,"%e\t",h[(int)i/100][(int)fmod(i,100)]);
			if(i<=*ntp)
			{
				fprintf(fa,"%e\t%e\n",g[(int)i/100][(int)fmod(i,100)],a[(int)i/100][(int)fmod(i,100)]);
			}
			else
			{
				fprintf(fa,"\n");
			}
		}
		fclose(fa);
	}
	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	for(i=0;i<3;i++)
	{
		free(xnu[i]);
	}
	free(xnu);
	free(rn);
	free(u);
	free(x);
	free(anyu);
	free(dife);

	return 0;
}

int nonline(int nelem,int ie,int *ndim,double **a,int **edge,int *nedge,int *nedges,double *anyu,double *dife,
		    int ngauss,double *xi,int **nod,double **xyz,int npoint,int *nvert,double **aa)
{
	int i,j,k,je,nod1,nod2,ned,iedges;
	double bhx,bh,fj,bxx,byy,bzz,bbb;
	double pi,df,dfx;
	double *u,*x,**dx,**dxi,**xn,**rn;

	iedges=*nedges;

	pi=acos(-1.0);

	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{       
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}
	
	for(i=0;i<ngauss;i++)
	{
		u[0]=xi[i];
		for(j=0;j<ngauss;j++)
		{
			u[1]=xi[j];
			for(k=0;k<ngauss;k++)
			{
				u[2]=xi[k];
				fj=fjacob(0,22,ie-1,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
				bxx=0.0;
				byy=0.0;
				bzz=0.0;
				for(je=0;je<*nedges;je++)
				{
					nod1=nod[inode11[je][0]][ie-1];
					nod2=nod[inode11[je][1]][ie-1];
					if(nod1>nod2)
					{
						ned=edge[ie-1][je]+1;
						aa[(int)ned/100][(int)fmod(ned,100)]=-a[(int)ned/100][(int)fmod(ned,100)];
					}
					if(nod2>nod1)
					{
						ned=edge[ie-1][je]+1;
						aa[(int)ned/100][(int)fmod(ned,100)]=a[(int)ned/100][(int)fmod(ned,100)];
					}
					bxx=bxx+rn[je][0]*aa[(int)ned/100][(int)fmod(ned,100)];
					byy=byy+rn[je][1]*aa[(int)ned/100][(int)fmod(ned,100)];
					bzz=bzz+rn[je][2]*aa[(int)ned/100][(int)fmod(ned,100)];
				}

				bbb=bxx*bxx+byy*byy+bzz*bzz;
				if(bbb<1e-20)
				{
					bbb=1e-20;
				}
				ss400(bbb,&bh,&df);

				bhx=bh;
				dfx=fabs(df);

				anyu[4*i+2*j+k]=bhx;
				dife[4*i+2*j+k]=dfx;
			}
		}
	}

	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	free(rn);
	free(u);
	free(x);

	return 0;
}

int banda(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,
		  int *nedges,int **ifa,double **a,int *nedge,double **fdb,double *bo,int *ntp,double **g,double **aa,
		  int jihn,double *xi,int **ifaf,int *nusef)
{
	int i,im,ned,ne,ie,je,nod1,nod2,iedges;
	double alpha,fj,fdbx,fdby,fdbz;
	double *u,*x,**dx,**dxi,**xn,**rn;
	FILE *fa,*fb,*fc;

	iedges=*nedges;

	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}

	alpha=1.0;
	fa=fopen("banda.dat","w");
	for(ned=1;ned<=*nedge;ned++)
	{
		im=ifa[(int)ned/100][(int)fmod(ned,100)];
		if(im>0)
		{
			a[(int)ned/100][(int)fmod(ned,100)]=a[(int)ned/100][(int)fmod(ned,100)]+g[(int)im/100][(int)fmod(im,100)]*alpha;
			fprintf(fa,"%d\t%d\t%e\t%e\n",ned,im,g[(int)im/100][(int)fmod(im,100)],a[(int)ned/100][(int)fmod(ned,100)]);
		}
	}
	fclose(fa);

	for(ne=0;ne<nelem;ne++)
	{
		bo[ne]=sqrt(fdb[ne][0]*fdb[ne][0]+fdb[ne][1]*fdb[ne][1]+fdb[ne][2]*fdb[ne][2]);
	}
	u[0]=0.0;
	u[1]=0.0;
	u[2]=0.0;
	fc=fopen("bandarn.dat","w");
	for(ie=0;ie<nelem;ie++)
	{
		fj=fjacob(0,22,ie,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
		fdbx=0.0;
		fdby=0.0;
		fdbz=0.0;
		fprintf(fc,"nelem(%d)\n",ie);
		fprintf(fc,"%e\t%e\t%e\n",rn[0][0],rn[0][1],rn[0][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[1][0],rn[1][1],rn[1][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[2][0],rn[2][1],rn[2][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[3][0],rn[3][1],rn[3][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[4][0],rn[4][1],rn[4][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[5][0],rn[5][1],rn[5][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[6][0],rn[6][1],rn[6][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[7][0],rn[7][1],rn[7][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[8][0],rn[8][1],rn[8][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[9][0],rn[9][1],rn[9][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[10][0],rn[10][1],rn[10][2]);
		fprintf(fc,"%e\t%e\t%e\n",rn[11][0],rn[11][1],rn[11][2]);
		for(je=0;je<*nedges;je++)
		{
			nod1=nod[inode11[je][0]][ie];
			nod2=nod[inode11[je][1]][ie];
			if(nod1>nod2)
			{
				ned=edge[ie][je]+1;
				aa[(int)ned/100][(int)fmod(ned,100)]=-a[(int)ned/100][(int)fmod(ned,100)];
			}
			if(nod2>nod1)
			{
				ned=edge[ie][je]+1;
				aa[(int)ned/100][(int)fmod(ned,100)]=a[(int)ned/100][(int)fmod(ned,100)];
			}
			fdbx=fdbx+rn[je][0]*aa[(int)ned/100][(int)fmod(ned,100)];
			fdby=fdby+rn[je][1]*aa[(int)ned/100][(int)fmod(ned,100)];
			fdbz=fdbz+rn[je][2]*aa[(int)ned/100][(int)fmod(ned,100)];
		}
		fdb[ie][0]=fdbx;
		fdb[ie][1]=fdby;
		fdb[ie][2]=fdbz;
	}
	fclose(fc);
	fb=fopen("bandafdb.dat","w");
	for(i=0;i<nelem;i++)
	{
		fprintf(fb,"%d\t%e\t%e\t%e\n",i,fdb[i][0],fdb[i][1],fdb[i][2]);
	}
	fclose(fb);

	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	free(rn);
	free(u);
	free(x);

	return 0;
}

int conver(int nelem,int *ndim,double *dmin,double **fdb,double *bo,int nficcg,int iihn,int *ifcon,FILE *ficcg)
{
	int ne,ifc;
	double wnr,bold,bnew,bdelta;

	ifc=0;
	wnr=0.0;

	for(ne=0;ne<nelem;ne++)
	{
		bold=bo[ne];
		bnew=sqrt(fdb[ne][0]*fdb[ne][0]+fdb[ne][1]*fdb[ne][1]+fdb[ne][2]*fdb[ne][2]);
		bdelta=fabs(bnew-bold);
		if(bdelta>*dmin)
		{
			ifc=ifc+1;
			wnr=wnr+bdelta;
		}
	}

	printf("iteration:%d\titako:%d\terrsum:%e\n",iihn,ifc,wnr);
/*	fprintf(ficcg,"iteration:%d\titako:%d\terrsum:%e\n",iihn,ifc,wnr);*/

	return 0;
}

int eddy(int **nod,int nelem,int *nvert,double **xyz,int npoint,int *ndim,int **edge,int *nedges,
		 int **imate,int *nblok,double **amate,int *nmate,double *freq,double **a,double **ao,int *nedge,int nfeddy,
		 double *dt,int *nusef,FILE *feddy)
{
	int ie,je,ibl,nod1,nod2,ned,i,iedges;
	double delta,cond,ejx,ejy,ejz,fj,efx,efy,efz;
	double *u,*x,**dx,**dxi,**xn,**rn,**dnf,aaa,aao;

	iedges=*nedges;

	u=(double*)malloc(sizeof(double)*3);
	x=(double*)malloc(sizeof(double)*3);
	dx=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dx[i]=(double*)malloc(sizeof(double)*3);
	}
	dxi=(double**)malloc(sizeof(double*)*3);
	for(i=0;i<3;i++)
	{
		dxi[i]=(double*)malloc(sizeof(double)*3);
	}
	xn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		xn[i]=(double*)malloc(sizeof(double)*3);
	}
	rn=(double**)malloc(sizeof(double*)*iedges);
	for(i=0;i<iedges;i++)
	{
		rn[i]=(double*)malloc(sizeof(double)*3);
	}
	dnf=(double**)malloc(sizeof(double*)*8);
	for(i=0;i<8;i++)
	{
		dnf[i]=(double*)malloc(sizeof(double)*3);
	}

	u[0]=0.0;
	u[1]=0.0;
	u[2]=0.0;

	for(ie=1;ie<=nelem;ie++)
	{
		for(ibl=0;ibl<*nblok;ibl++)
		{
			if((ie>=imate[ibl][0])&&(ie<=imate[ibl][1]))
			{
				delta=amate[imate[ibl][2]-1][3]/(*dt);
				cond=amate[imate[ibl][2]-1][3];
			}
		}
		fj=fjacob(0,14,ie-1,u,x,dx,dxi,xn,rn,nod,xyz,npoint,nelem,nvert,ndim,nedges);
		ejx=0.0;
		ejy=0.0;
		ejz=0.0;
		for(je=0;je<*nedges;je++)
		{
			nod1=nod[inode11[je][0]][ie-1];
			nod2=nod[inode11[je][1]][ie-1];
			ned=edge[ie-1][je]+1;
			if(nod1>nod2)
			{
				aao=-ao[(int)ned/100][(int)fmod(ned,100)];
				aaa=-a[(int)ned/100][(int)fmod(ned,100)];
			}
			if(nod2>nod1)
			{
				aao=ao[(int)ned/100][(int)fmod(ned,100)];
				aaa=a[(int)ned/100][(int)fmod(ned,100)];
			}
			ejx=ejx-delta*xn[je][0]*(aaa-aao);
			ejy=ejy-delta*xn[je][1]*(aaa-aao);
			ejz=ejz-delta*xn[je][2]*(aaa-aao);
		}
		if(*nusef==1)
		{
			efx=0.0;
			efy=0.0;
			efz=0.0;
		}
		fprintf(feddy,"%e\t%e\t%e\n",ejx,ejy,ejz);
	}

	for(i=0;i<3;i++)
	{
		free(dx[i]);
	}
	free(dx);
	for(i=0;i<3;i++)
	{
		free(dxi[i]);
	}
	free(dxi);
	for(i=0;i<iedges;i++)
	{
		free(xn[i]);
	}
	free(xn);
	for(i=0;i<iedges;i++)
	{
		free(rn[i]);
	}
	free(rn);
	for(i=0;i<8;i++)
	{
		free(dnf[i]);
	}
	free(dnf);
	free(u);
	free(x);

	
	return 0;
}

int atoao(double **a,double **ao,int *nedge)
{
	int ned;

	for(ned=1;ned<=*nedge;ned++)
	{
		ao[(int)ned/100][(int)fmod(ned,100)]=a[(int)ned/100][(int)fmod(ned,100)];
	}

	return 0;
}

/* ICCG Method Start */

int iccg(int *mband,int *mun,double **h,double **g,double **al,int **nadj,int **nsky,int **nfront,int *nitecg,double *epsicg,
		 int nitenr,int istep,int nficcg,FILE *ficcg,int npoint,double **s,double **a,double **r,double **p)
{
/*	���v�ϐ�
	nband	:��0�_�̐�
	nun		:�W���s���̃T�C�Y
	nfront	:��0�_�̏W��(size:nband)
	nadj	:�Ƃ����ӗv�f�ɑ��݂�����0�_�̐�(size:nun)
	nsky	:�Ƃ����ӗv�f�̔�0�_��nfront�̂ǂ̏ꏊ�����n�܂邩���L�^(size:nun)
	r		:�c���x�N�g��
	g		:���Ӄx�N�g���y�ѕ⏕�x�N�g��
	h		:�W���s��
	s		:�X�P�[�����O�t�@�N�^
	p		:�w�������x�N�g��
	a		:���x�N�g��																	*/
	int nband,nun,iscale,iconv1,iconv2;
	int i,l,nic;
	double sign,aegmax,gmax,g2,cnew,r2;
	double bet,cold,alf,bunbo,relerr,eps;
	double **q;
	FILE *fa,*fb,*fc;
	double t1,t2,t3,t4,t5,tic,tcg;
	t1=gettime();
	fc=fopen("timeiccg.dat","w");
	nband=*mband;
	nun=*mun;
	fa=fopen("nfront.dat","w");
	for(i=1;i<=nband;i++)
	{
		fprintf(fa,"%d\t",i);
		if(i<=nun)
		{
			fprintf(fa,"%d\t%d\t",nsky[(int)i/100][(int)fmod(i,100)],nadj[(int)i/100][(int)fmod(i,100)]);
		}
		else
		{
			fprintf(fa,"\t\t");
		}
		fprintf(fa,"%d\n",nfront[(int)i/100][(int)fmod(i,100)]);
	}
	fclose(fa);
	double t0=gettime();
	if(nun<=0)
	{
		printf("iccg errer! --- nun <= 0 ---\n");
		nic=0;
		return 1;
	}

	printf("Now be in subroutine-ICCG method!\n");

	iscale=1;

	symsca(nband,nun,h,a,g,s,nsky,nadj,nfront,iscale);

	inchde(nband,nun,h,al,nadj,nsky,nfront,nficcg,ficcg);
/*�����@*/
/*	rep_inchde(nband,nun,h,al,nadj,nsky,nfront);*/
fa=fopen("iccgv1.dat","w");
	forbac(nband,nun,al,g,a,nadj,nsky,nfront);
	fb=fopen("icend.dat","w");
	for(i=1;i<=nband;i++)
	{
		fprintf(fb,"%e\n",al[(int)i/100][(int)fmod(i,100)]);
	}
	fclose(fb);
	for(i=1;i<=nun;i++)
	{
		r[(int)i/100][(int)fmod(i,100)]=g[(int)i/100][(int)fmod(i,100)];
	}

	sign=-1.0;

	multi(nband,nun,h,a,r,nadj,nsky,nfront,sign);

	maxval(nun,g,&gmax);
	eps=*epsicg;
	aegmax=eps*gmax;
	printf("aegmax= %e\n",aegmax);
	produc(nun,g,g,&g2);
	g2=sqrt(g2);
	printf("g2= %e\n",g2);

	q=(double**)malloc(sizeof(double*)*(int)(npoint/25+1));
	for(i=0;i<npoint/25+1;i++)
	{
		q[i]=(double*)malloc(sizeof(double)*100);
	}

	t2=gettime();
	tic=0.0;
	tcg=0.0;
	for(l=1;l<=*nitecg;l++)
	{
		t3=gettime();
		forbac(nband,nun,al,r,g,nadj,nsky,nfront);
		t4=gettime();
		produc(nun,r,g,&cnew);

		if(l<2)
		{
			for(i=1;i<=nun;i++)
			{
				p[(int)i/100][(int)fmod(i,100)]=g[(int)i/100][(int)fmod(i,100)];
			}
		}
		else
		{
			bet=cnew/cold;
			for(i=1;i<=nun;i++)
			{
				p[(int)i/100][(int)fmod(i,100)]=g[(int)i/100][(int)fmod(i,100)]+bet*p[(int)i/100][(int)fmod(i,100)];
			}
		}
		cold=cnew;
		for(i=1;i<=nun;i++)
		{
			q[(int)i/100][(int)fmod(i,100)]=0.0;
		}

		sign=1.0;

		multi(nband,nun,h,p,q,nadj,nsky,nfront,sign);

		produc(nun,p,q,&bunbo);

		iconv1=0;
		iconv2=0;

		if(fabs(bunbo)<1e-30)
		{
			printf("fabs(bunbo) = %e <1.0e-30\n",fabs(bunbo));
			if(nficcg>0)
			{
/*				fprintf(ficcg,"fabs(bunbo) = %e <1.0e-30\n",fabs(bunbo));*/
			}
			goto lend;
		}
		alf=cnew/bunbo;
		for(i=1;i<=nun;i++)
		{
			a[(int)i/100][(int)fmod(i,100)]=a[(int)i/100][(int)fmod(i,100)]+alf*p[(int)i/100][(int)fmod(i,100)];
			r[(int)i/100][(int)fmod(i,100)]=r[(int)i/100][(int)fmod(i,100)]-alf*q[(int)i/100][(int)fmod(i,100)];
			if(fabs(r[(int)i/100][(int)fmod(i,100)])>aegmax)
			{
				iconv1=iconv1+1;
			}
		}

		produc(nun,r,r,&r2);
		r2=sqrt(r2);

		relerr=r2/g2;
		fprintf(fa,"%d\t%d\n",l,iconv1);
		if(relerr>*epsicg)
		{
			iconv2=1;
		}

		if(fmod(l,20)==0)
		{
			printf(" l=%d\n iconv1=%d\n epsicg=%e\n relerr=%e\n",l,iconv1,*epsicg,relerr);
		}
		t5=gettime();
		tic+=(t4-t3);
		tcg+=(t5-t4);
		if(iconv1==0)
		{
			goto lend;
		}
	}
lend:;
fclose(fa);
	iscale=2;
	tic+=(t2-t1);
	fprintf(fc,"%e\t%e\n",tic,tcg);
	fclose(fc);
	symsca(nband,nun,h,a,q,s,nsky,nadj,nfront,iscale);

	for(i=1;i<=nun;i++)
	{
		g[(int)i/100][(int)fmod(i,100)]=a[(int)i/100][(int)fmod(i,100)];
	}

	printf("the number of iterations = %d\n",l);
	if(nficcg>0)
	{
/*		fprintf(ficcg,"the number of iterations = %d\n",l);*/
	}
	nic=l;

	if((iconv1==0)&&(iconv2==0))
	{
		printf("----- convergence OK ! -----\n");
		if(nficcg>0)
		{
/*			fprintf(ficcg,"----- convergence OK ! -----\n");*/
		}
	}
	else
	{
		printf("----- not convergence ! -----\n");
		if(nficcg>0)
		{
/*			fprintf(ficcg,"----- not convergence ! -----\n");*/
		}
	}
	fb=fopen("iccgr.dat","w");
	for(i=1;i<=nun;i++)
	{
		fprintf(fb,"%e",g[(int)i/100][(int)fmod(i,100)]);
		fprintf(fb,"\n");
	}
	fclose(fb);
	for(i=0;i<(int)(npoint/25+1);i++)
	{
		free(q[i]);
	}
	free(q);
	

	return 0;
}

int symsca(int nband,int nun,double **h,double **a,double **g,double **s,int **nsky,int **nadj,int **nfront,int iscale)
{
	int i,j,ins,ina,ii,nij,nf;
	double hii;
	FILE *fa;

	if(iscale!=2)
	{
		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)];
			ii=ins+ina;
			hii=sqrt(fabs(h[(int)ii/100][(int)fmod(ii,100)]));
			s[(int)i/100][(int)fmod(i,100)]=1.0/hii;
			g[(int)i/100][(int)fmod(i,100)]=g[(int)i/100][(int)fmod(i,100)]*s[(int)i/100][(int)fmod(i,100)];
		}
		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)];
			for(j=1;j<=ina;j++)
			{
				nij=ins+j;
				nf=nfront[(int)nij/100][(int)fmod(nij,100)];
				h[(int)nij/100][(int)fmod(nij,100)]=h[(int)nij/100][(int)fmod(nij,100)]*s[(int)i/100][(int)fmod(i,100)]*s[(int)nf/100][(int)fmod(nf,100)];
			}
		}
	}
	else
	{
		for(i=1;i<=nun;i++)
		{
			a[(int)i/100][(int)fmod(i,100)]=a[(int)i/100][(int)fmod(i,100)]*s[(int)i/100][(int)fmod(i,100)];
		}
	}
	fa=fopen("h.dat","w");
	for(i=1;i<=nun;i++)
	{
		fprintf(fa,"%d\t%e\t%e\t%e\t%e\n",i,h[(int)i/100][(int)fmod(i,100)],s[(int)i/100][(int)fmod(i,100)],a[(int)i/100][(int)fmod(i,100)],g[(int)i/100][(int)fmod(i,100)]);
	}
	fclose(fa);

	return 0;
}

int inchde(int nband,int nun,double **h,double **al,int **nadj,int **nsky,int **nfront,int nficcg,FILE *ficcg)
{
	int icount,i,ins,ina,nid,j,nij;
	int jns,jnf,jna,kmx,nik,knf,knd;
	int nil,k,l,lnf,jnd;
	double alpha;
	FILE *fa;

	alpha=1.0;
	for(icount=0;icount<100;icount++)
	{
		printf("alpha[%d] = %e\n",icount,alpha);
		if(nficcg>0)
		{
/*			fprintf(ficcg,"alpha[%d] = %e\n",icount,alpha);*/
		}
		for(i=1;i<=nband;i++)
		{
			al[(int)i/100][(int)fmod(i,100)]=h[(int)i/100][(int)fmod(i,100)];
		}
		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)];
			nid=ins+ina;
/*			printf("%d\t%d\t%d\t%d\t%e\n",i,ins,ina,nid,h[(int)nid/100][(int)fmod(nid,100)]);	*/
			al[(int)nid/100][(int)fmod(nid,100)]=h[(int)nid/100][(int)fmod(nid,100)]*alpha;
		}
		fa=fopen("inchde.dat","w");
		for(i=1;i<=nband;i++)
		{
			fprintf(fa,"%d\t%e\n",i,al[(int)i/100][(int)fmod(i,100)]);
		}
		fclose(fa);
		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)]-1;
			nid=ins+ina+1;
			for(j=2;j<=ina;j++)
			{
				nij=ins+j;
				jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
				jns=nsky[(int)jnf/100][(int)fmod(jnf,100)]-1;
				jna=nadj[(int)jnf/100][(int)fmod(jnf,100)]-1;
				kmx=j-1;
				for(k=1;k<=kmx;k++)
				{
					nik=ins+k;
					knf=nfront[(int)nik/100][(int)fmod(nik,100)];
					knd=nsky[(int)knf/100][(int)fmod(knf,100)]+nadj[(int)knf/100][(int)fmod(knf,100)]-1;
					for(l=1;l<=jna;l++)
					{
						nil=jns+l;
						lnf=nfront[(int)nil/100][(int)fmod(nil,100)];
						if(lnf==knf)
						{
							al[(int)nij/100][(int)fmod(nij,100)]=al[(int)nij/100][(int)fmod(nij,100)]-al[(int)nik/100][(int)fmod(nik,100)]*al[(int)nil/100][(int)fmod(nil,100)]*al[(int)knd/100][(int)fmod(knd,100)];
						}
					}
				}
			}
			for(j=1;j<=ina;j++)
			{
				nij=ins+j;
				jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
				jnd=nsky[(int)jnf/100][(int)fmod(jnf,100)]+nadj[(int)jnf/100][(int)fmod(jnf,100)]-1;
				al[(int)nid/100][(int)fmod(nid,100)]=al[(int)nid/100][(int)fmod(nid,100)]-al[(int)nij/100][(int)fmod(nij,100)]*al[(int)nij/100][(int)fmod(nij,100)]*al[(int)jnd/100][(int)fmod(jnd,100)];
			}
			al[(int)nid/100][(int)fmod(nid,100)]=1.0/al[(int)nid/100][(int)fmod(nid,100)];
			if(al[(int)nid/100][(int)fmod(nid,100)]<0.0)
			{
				goto Warn;
			}
		}
		goto Clear;
Warn:;
		printf("Warning d(%d) = %e\n",i,al[(int)nid/100][(int)fmod(nid,100)]);
		if(nficcg>0)
		{
	/*		fprintf(ficcg,"Warning d(%d) = %e\n",nid,al[(int)nid/100][(int)fmod(nid,100)]);*/
		}
		alpha=alpha+0.01;
		if(icount==99)
		{
			alpha=2.0;
		}
	}
Clear:;
	printf("icount = %d\n",icount);
	if(nficcg>0)
	{
/*		fprintf(ficcg,"icount = %d\n",icount);*/
	}
	return 0;
}

int forbac(int nband,int nun,double **al,double **r,double **q,int **nadj,int **nsky,int **nfront)
{

	forsub(nband,nun,al,r,q,nadj,nsky,nfront);

	bacsub(nband,nun,al,q,nadj,nsky,nfront);

	return 0;
}

int forsub(int nband,int nun,double **al,double **r,double **v,int **nadj,int **nsky,int **nfront)
{
	int i,j,ins,ina,nij,nnf,jnf,jnd;
	double ali;

	for(i=1;i<=nun;i++)
	{
		ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
		ina=nadj[(int)i/100][(int)fmod(i,100)]-1;
		ali=al[(int)(ins+ina+1)/100][(int)fmod(ins+ina+1,100)];
		v[(int)i/100][(int)fmod(i,100)]=r[(int)i/100][(int)fmod(i,100)]*ali;
	}
	for(i=2;i<=nun;i++)
	{
		ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
		ina=nadj[(int)i/100][(int)fmod(i,100)]-1;
		ali=al[(int)(ins+ina+1)/100][(int)fmod(ins+ina+1,100)];
		for(j=1;j<=ina;j++)
		{
			nij=ins+j;
			nnf=nfront[(int)nij/100][(int)fmod(nij,100)];
			v[(int)i/100][(int)fmod(i,100)]=v[(int)i/100][(int)fmod(i,100)]-al[(int)nij/100][(int)fmod(nij,100)]*v[(int)nnf/100][(int)fmod(nnf,100)]*ali;
		}
	}
	return 0;
}
int bacsub(int nband,int nun,double **al,double **q,int **nadj,int **nsky,int **nfront)
{
	int i,j,ins,ina,nij,nnf,jnf,jnd;
	double ali;

	for(i=nun;i>=2;)
	{
		ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
		ina=nadj[(int)i/100][(int)fmod(i,100)]-1;
		for(j=1;j<=ina;j++)
		{
			nij=ins+j;
			jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
			jnd=nsky[(int)jnf/100][(int)fmod(jnf,100)]+nadj[(int)jnf/100][(int)fmod(jnf,100)]-1;
			q[(int)jnf/100][(int)fmod(jnf,100)]=q[(int)jnf/100][(int)fmod(jnf,100)]-al[(int)nij/100][(int)fmod(nij,100)]*al[(int)jnd/100][(int)fmod(jnd,100)]*q[(int)i/100][(int)fmod(i,100)];
		}
		i=i-1;
	}
	return 0;
}
int produc(int nun,double **r,double **q,double *cnew)
{
	int i;
	double cn;
	cn=0.0e0;
	for(i=1;i<=nun;i++)
	{
		cn=cn+r[(int)i/100][(int)fmod(i,100)]*q[(int)i/100][(int)fmod(i,100)];
	}
	*cnew=cn;
	return 0;
}

int multi(int nband,int nun,double **h,double **p,double **g,int **nadj,int **nsky,int **nfront,double sign)
{
	int i,j,ins,ina,nij,jp,jnf;
	double pi;

	for(i=1;i<=nun;i++)
	{
		ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
		ina=nadj[(int)i/100][(int)fmod(i,100)]-1;
		for(j=1;j<=ina;j++)
		{
			nij=ins+j;
			jp=nfront[(int)nij/100][(int)fmod(nij,100)];
			g[(int)i/100][(int)fmod(i,100)]=g[(int)i/100][(int)fmod(i,100)]+sign*h[(int)nij/100][(int)fmod(nij,100)]*p[(int)jp/100][(int)fmod(jp,100)];
		}
	}
	for(i=1;i<=nun;i++)
	{
		ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
		ina=nadj[(int)i/100][(int)fmod(i,100)];
		pi=p[(int)i/100][(int)fmod(i,100)];
		for(j=1;j<=ina;j++)
		{
			nij=ins+j;
			jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
			g[(int)jnf/100][(int)fmod(jnf,100)]=g[(int)jnf/100][(int)fmod(jnf,100)]+sign*pi*h[(int)nij/100][(int)fmod(nij,100)];
		}
	}

	return 0;
}

int maxval(int nun,double **g,double *gmax)
{
	int i;
	double gm;

	gm=0.0;
	for(i=1;i<=nun;i++)
	{
		if(fabs(g[(int)i/100][(int)fmod(i,100)])>gm)
		{
			gm=fabs(g[(int)i/100][(int)fmod(i,100)]);
		}
	}
	*gmax=gm;

	return 0;
}

int rep_inchde(int nband,int nun,double **h,double **al,int **nadj,int **nsky,int **nfront)
{
	int icount,i,ins,ina,nid,j,nij,icon;
	int jns,jnf,jna,nik,knf,knd,jnd;
	int nil,k,l,lnf;
	double eps,omg;
	double **alo;
	double **r;
	FILE *fa,*fb;
	eps=0.09;
	omg=1.8;
	alo=(double**)malloc(sizeof(double*)*((int)nband/100+1));
	for(i=0;i<nband/100+1;i++)
	{
		alo[i]=(double*)malloc(sizeof(double)*100);
	}
	r=(double**)malloc(sizeof(double*)*((int)nband/100+1));
	for(i=0;i<nband/100+1;i++)
	{
		r[i]=(double*)malloc(sizeof(double)*100);
	}

	for(i=1;i<=nband;i++)
	{
		alo[(int)i/100][(int)fmod(i,100)]=0.0;
	}
	for(i=1;i<=nun;i++)
	{
		nid=nsky[(int)i/100][(int)fmod(i,100)]+nadj[(int)i/100][(int)fmod(i,100)]-1;
		alo[(int)nid/100][(int)fmod(nid,100)]=sqrt(fabs(h[(int)nid/100][(int)fmod(nid,100)]));
	}
	icount=0;
	fa=fopen("conv.dat","w");
	
	while(1)
	{
		icount++;

		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)];
			nid=ins+ina;
			for(j=1;j<=ina;j++)
			{
				nij=ins+j;
				jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
				jnd=nsky[(int)jnf/100][(int)fmod(jnf,100)]+nadj[(int)jnf/100][(int)fmod(jnf,100)]-1;
				if(nij==nid)
				{
					al[(int)nij/100][(int)fmod(nij,100)]=h[(int)nij/100][(int)fmod(nij,100)];
				}
				else
				{
					al[(int)nij/100][(int)fmod(nij,100)]=h[(int)nij/100][(int)fmod(nij,100)]/alo[(int)jnd/100][(int)fmod(jnd,100)];
				}
				r[(int)nij/100][(int)fmod(nij,100)]=h[(int)nij/100][(int)fmod(nij,100)];
			}

		}
		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)]-1;
			nid=ins+ina+1;
			for(j=1;j<=ina;j++)
			{
				nij=ins+j;
				jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
				jns=nsky[(int)jnf/100][(int)fmod(jnf,100)]-1;
				jna=nadj[(int)jnf/100][(int)fmod(jnf,100)]-1;
				jnd=jns+jna+1;
				for(k=1;k<=j-1;k++)
				{
					nik=ins+k;
					knf=nfront[(int)nik/100][(int)fmod(nik,100)];
					knd=nsky[(int)knf/100][(int)fmod(knf,100)]+nadj[(int)knf/100][(int)fmod(knf,100)]-1;
					for(l=1;l<=jna;l++)
					{
						nil=jns+l;
						lnf=nfront[(int)nil/100][(int)fmod(nil,100)];
						if(lnf==knf)
						{
							al[(int)nij/100][(int)fmod(nij,100)]=al[(int)nij/100][(int)fmod(nij,100)]-alo[(int)nik/100][(int)fmod(nik,100)]*alo[(int)nil/100][(int)fmod(nil,100)]/alo[(int)jnd/100][(int)fmod(jnd,100)];
						}
					}
				}
			}
			for(j=1;j<=ina;j++)
			{
				nij=ins+j;
				al[(int)nid/100][(int)fmod(nid,100)]=al[(int)nid/100][(int)fmod(nid,100)]-alo[(int)nij/100][(int)fmod(nij,100)]*alo[(int)nij/100][(int)fmod(nij,100)];
			}
			al[(int)nid/100][(int)fmod(nid,100)]=sqrt(fabs(al[(int)nid/100][(int)fmod(nid,100)]));
		}
		icon=0;
		for(i=1;i<=nun;i++)
		{
			ins=nsky[(int)i/100][(int)fmod(i,100)]-1;
			ina=nadj[(int)i/100][(int)fmod(i,100)];
			nid=ins+ina;
			for(j=1;j<=ina;j++)
			{
				nij=ins+j;
				jnf=nfront[(int)nij/100][(int)fmod(nij,100)];
				jns=nsky[(int)jnf/100][(int)fmod(jnf,100)]-1;
				jna=nadj[(int)jnf/100][(int)fmod(jnf,100)];
				for(k=1;k<=j;k++)
				{
					nik=ins+k;
					knf=nfront[(int)nik/100][(int)fmod(nik,100)];
					knd=nsky[(int)knf/100][(int)fmod(knf,100)]+nadj[(int)knf/100][(int)fmod(knf,100)]-1;
					for(l=1;l<=jna;l++)
					{
						nil=jns+l;
						lnf=nfront[(int)nil/100][(int)fmod(nil,100)];
						if(lnf==knf)
						{
							r[(int)nij/100][(int)fmod(nij,100)]=r[(int)nij/100][(int)fmod(nij,100)]-al[(int)nik/100][(int)fmod(nik,100)]*al[(int)nil/100][(int)fmod(nil,100)];
						}
					}
				}
				if(fabs(r[(int)nij/100][(int)fmod(nij,100)])>eps)
				{
					icon=icon+1;
				}
			}
		}
		fprintf(fa,"%d\t%d\n",icount,icon);
		if(icon==0)
		{
			goto clear;
		}
		for(i=1;i<=nband;i++)
		{
			alo[(int)i/100][(int)fmod(i,100)]=al[(int)i/100][(int)fmod(i,100)];
		}
		if(icount==5)
		{
			printf("not convergence!\n");
			goto errered;
		}
	}
errered:;
	fb=fopen("conv_errer.dat","w");
	for(i=1;i<=nband;i++)
	{
		if(fabs(r[(int)i/100][(int)fmod(i,100)])>eps)
		{
			fprintf(fb,"%d\t%e\n",i,r[(int)i/100][(int)fmod(i,100)]);
		}
	}
	fclose(fb);
clear:;
	fclose(fa);
	for(i=1;i<=nun;i++)
	{
		nid=nsky[(int)i/100][(int)fmod(i,100)]+nadj[(int)i/100][(int)fmod(i,100)]-1;
		al[(int)nid/100][(int)fmod(nid,100)]=1.0/al[(int)nid/100][(int)fmod(nid,100)];
	}
	for(i=0;i<(int)nband/100+1;i++)
	{
		free(alo[i]);
	}
	free(alo);
	for(i=0;i<(int)nband/100+1;i++)
	{
		free(r[i]);
	}
	free(r);
	return 0;
}

int bhodf(double bb,double *bh,double *df)
{
	int init,i;
	double ams,anyu,b;
	double bw[38]={0.0000,
                   0.0025,
                   0.0050,
                   0.0125,
                   0.0250,
                   0.0500,
                   0.1000,
                   0.2000,
                   0.3000,
                   0.4000,
                   0.5000,
                   0.6000,
                   0.7000,
                   0.8000,
                   0.9000,      
                   1.0000,
                   1.1000,
                   1.2000,
                   1.3000,
                   1.4000,
                   1.5000,
                   1.5500,
                   1.6000,
                   1.6500,
                   1.7000,
                   1.7500,
                   1.8000,
                   1.8760,
                   1.9446,
                   2.0059,      
                   2.0598,
                   2.1063,
                   2.1454,
                   2.1772,
                   2.2016,
                   2.2186,
                   2.2282,
	               2.2342};

	double ch[38][4]={{2.02579E+08,-5.06448E+05, 7.23960E+03, 0.00000E+00,},
	{-5.51168E+08, 1.49467E+06, 8.50572E+03, 1.80990E+01},
	{3.56552E+07,-7.30544E+05, 5.64468E+03, 4.00930E+01},
	{6.39962E+06,-9.98773E+04, 7.03327E+02, 5.63770E+01},
	{-1.53354E+06, 5.31588E+04, 1.20621E+03, 6.20620E+01},
	{-7.67461E+04,-1.32996E+03, 9.88763E+02, 1.01480E+02},
	{-3.21997E+03,-1.05719E+02, 2.80172E+02, 1.38000E+02},
	{-2.09522E+03, 1.55233E+02, 1.62429E+02, 1.61740E+02},
	{-1.80096E+03, 8.69090E+01, 1.30619E+02, 1.77440E+02},
	{-1.04622E+01,-9.67102E+00, 9.39717E+01, 1.89570E+02},
	{4.69146E+01,-6.92825E+00, 9.17236E+01, 1.98860E+02},
	{-2.81992E+02, 7.97460E+01, 9.17454E+01, 2.08010E+02},
	{2.02468E+03, 6.81839E+01, 9.92348E+01, 2.17700E+02},
	{-1.89232E+03, 3.42112E+02, 1.73612E+02, 2.30330E+02},
	{6.23073E+03,-6.81723E+02, 1.85265E+02, 2.49220E+02},
	{-4.27411E+03, 8.72988E+02, 2.35842E+02, 2.67160E+02},
	{1.68362E+02, 3.99960E+01, 2.82216E+02, 2.95200E+02},
	{-2.37323E+04, 5.45957E+03, 2.95267E+02, 3.23990E+02},
	{-1.07043E+04, 3.49930E+03, 6.75212E+02, 3.84380E+02},
	{-5.10605E+04, 1.45106E+04, 1.05394E+03, 4.76190E+02},
	{-3.08687E+05, 5.03415E+04, 2.42425E+03, 6.75630E+02},
	{-5.17697E+05, 9.38157E+04, 5.14324E+03, 8.84110E+02},
	{-4.67137E+05, 1.33076E+05, 1.06421E+04, 1.31110E+03},
	{-9.96949E+05, 1.91564E+05, 2.04461E+04, 2.11750E+03},
	{3.19159E+06,-2.40758E+03, 3.21254E+04, 3.49410E+03},
	{-6.10222E+05, 1.03961E+05, 5.58216E+04, 5.49330E+03},
	{1.74379E+05,-3.50925E+03, 6.16410E+04, 8.46800E+03},
	{-3.80159E+05, 9.86974E+04, 6.41292E+04, 1.32090E+04},
	{-2.94331E+05, 1.00222E+05, 7.23035E+04, 1.79500E+04},
	{-5.89985E+05, 1.55855E+05, 8.12726E+04, 2.26910E+04},
	{-1.19338E+06, 2.49574E+05, 9.29317E+04, 2.74320E+04},
	{-2.94421E+06, 4.43172E+05, 1.08401E+05, 3.21730E+04},
	{-8.23834E+06, 8.76266E+05, 1.29554E+05, 3.69130E+04},
	{-3.24408E+07, 2.18549E+06, 1.60291E+05, 4.16540E+04},
	{-7.66223E+07, 5.41327E+06, 2.09001E+05, 4.63950E+04},
	{-3.62341E+07, 1.77685E+07, 3.26621E+05, 5.11360E+04},
	{-3.67787E+09, 4.41347E+07, 6.57756E+05, 5.58770E+04},
	{0.00000E+00, 0.00000E+00, 7.90161E+05, 6.06180E+04}};

	init=0;
	ams=0.21580e1;
	anyu=795774.7155;

	bb=sqrt(bb);
	if(init==1)
	{
		goto bhodf1;
	}
	init=1;
bhodf1:;
	for(i=1;i<38;i++)
	{
		if(bb<bw[i])
		{
			goto bhodf2;
		}
	}
	*bh=anyu*(1.0-ams/bb);
	bb=bb*bb;
	*df=0.5/bb*(anyu-(*bh));
	return 0;
bhodf2:;
	i=i-1;
	b=bb-bw[i];
	*bh=((ch[i][0]*b+ch[i][1])*b+ch[i][2])*b+ch[i][3];
	*df=(3.0*ch[i][0]*b+2.0*ch[i][1])*b+ch[i][2];
	*bh=(*bh)/bb;
	bb=bb*bb;
	*df=0.5/bb*((*df)-(*bh));

	return 0;
}

int bhodf21(double bb,double *bh,double *df)
{
	int init,i;
	double ams,anyu,b;
	double bw[54]={0.0100,
                   0.0200,
                   0.0300,
                   0.0400,
                   0.0500,
                   0.0600,
                   0.0700,
                   0.0800,
                   0.0900,
                   0.1000,
                   0.1100,
                   0.1500,
                   0.2000,
                   0.3000,
                   0.4000,     
                   0.5000,
                   0.6000,
                   0.7000,
                   0.8000,
                   0.9000,
                   1.0000,
                   1.1000,
                   1.2000,
                   1.2500,
                   1.3000,
                   1.3200,
                   1.3400,
                   1.3600,
                   1.3800,
                   1.4000,    
                   1.4200,
                   1.4400,
                   1.4600,
                   1.4800,
                   1.5000,
                   1.5200,
                   1.5400,
                   1.5600,
                   1.5800,
                   1.6000,
                   1.6200,
                   1.6400,
                   1.6600,
                   1.6800,
                   1.7000,      
                   1.8000,
                   1.8200,
                   1.8600,
                   1.9878,
                   2.0898,
                   2.1659,
                   2.2162,
                   2.2406,
				   2.2521};

	double ch[54][4]={{ 1.66677E+05,-6.66677E+03, 1.75000E+03, 2.00000E+01},
	{-3.33322E+05,-3.33356E+03, 1.66667E+03, 3.70000E+01},
	{ 6.00000E+05,-1.60000E+04, 1.50000E+03, 5.30000E+01},
	{-8.99996E+05, 2.99995E+03, 1.36000E+03, 6.70000E+01},
	{ 7.49996E+05,-2.24999E+04, 1.15000E+03, 8.00000E+01},
	{-2.49997E+05,-3.05176E-02, 9.25000E+02, 9.00000E+01},
	{-4.27246E+00,-4.99995E+03, 8.50000E+02, 9.90000E+01},
	{-4.99992E+05,-9.76562E-02, 7.50000E+02, 1.07000E+02},
	{ 0.00000E+00, 0.00000E+00, 6.00000E+02, 1.14000E+02},
	{ 0.00000E+00, 0.00000E+00, 6.00000E+02, 1.20000E+02},
	{ 3.96636E+04,-5.33655E+03, 6.00000E+02, 1.26000E+02},
	{ 5.28846E+04,-3.91346E+03, 3.63461E+02, 1.44000E+02},
	{ 4.87501E+03,-7.50006E+01, 3.68750E+02, 1.59000E+02},
	{ 0.00000E+00, 0.00000E+00, 5.00000E+02, 2.00000E+02},
	{ 0.00000E+00, 0.00000E+00, 5.00000E+02, 2.50000E+02},
	{-1.49999E+04, 3.49999E+03, 5.00000E+02, 3.00000E+02},
	{ 4.99995E+03, 6.71387E-03, 7.50000E+02, 3.70000E+02},
	{-5.99995E+03, 1.59999E+03, 9.00000E+02, 4.50000E+02},
	{-1.71431E+03, 7.71435E+02, 1.04000E+03, 5.50000E+02},
	{-1.18681E+04, 3.75823E+03, 1.14286E+03, 6.60000E+02},
	{-8.65398E+03, 5.48078E+03, 1.53846E+03, 8.00000E+02},
	{ 3.75003E+04, 2.49997E+03, 2.37500E+03, 1.00000E+03},
	{ 0.00000E+00, 0.00000E+00, 4.00000E+03, 1.30000E+03},
	{ 0.00000E+00, 0.00000E+00, 4.00000E+03, 1.50000E+03},
	{-8.33218E+05, 6.66628E+04, 4.00000E+03, 1.70000E+03},
	{-2.08366E+05, 2.08348E+04, 5.66666E+03, 1.80000E+03},
	{ 1.87501E+06,-2.50004E+04, 6.25001E+03, 1.92000E+03},
	{ 0.00000E+00, 0.00000E+00, 7.50001E+03, 2.05000E+03},
	{ 0.00000E+00, 0.00000E+00, 7.50001E+03, 2.20000E+03},
	{ 0.00000E+00, 0.00000E+00, 7.50001E+03, 2.35000E+03}, 
	{-2.08318E+06, 1.66660E+05, 7.50001E+03, 2.50000E+03},
	{-1.45836E+07, 4.58342E+05, 1.16666E+04, 2.70000E+03},
	{ 1.87500E+07,-4.99999E+05, 1.25000E+04, 3.00000E+03},
	{ 0.00000E+00, 0.00000E+00, 1.50000E+04, 3.20000E+03},
	{ 0.00000E+00, 0.00000E+00, 1.50000E+04, 3.50000E+03},
	{ 0.00000E+00, 0.00000E+00, 1.50000E+04, 3.80000E+03},
	{-6.25030E+06, 3.75006E+05, 1.50000E+04, 4.10000E+03},
	{ 2.08357E+06, 8.33260E+04, 2.24999E+04, 4.50000E+03},
	{-1.04183E+06, 1.04174E+05, 2.83333E+04, 5.00000E+03},
	{ 3.12505E+06,-1.75781E+00, 3.12500E+04, 5.60000E+03},
	{-4.16668E+06, 2.08334E+05, 3.50000E+04, 6.25000E+03},
	{ 8.33478E+05, 6.66639E+04, 3.83334E+04, 7.00000E+03},
	{ 4.99837E+06, 3.00031E+05, 4.20001E+04, 7.80000E+03},
	{-4.49214E+02, 2.69530E+01, 5.99993E+04, 8.80000E+03},
	{ 5.39064E+01,-3.59375E+00, 5.99999E+04, 1.00000E+04},
	{-3.25315E+07, 1.15057E+06, 6.00008E+04, 1.60000E+04},
	{ 2.13239E+06,-1.34939E+05, 6.69858E+04, 1.74000E+04},
	{-1.19425E+05, 5.56560E+04, 6.64261E+04, 2.00000E+04},
	{-7.64854E+05, 2.23959E+05, 7.48002E+04, 2.91490E+04},
	{-2.62239E+06, 5.09789E+05, 9.66151E+04, 3.82970E+04},
	{-1.12853E+07, 1.62617E+06, 1.28644E+05, 4.74460E+04},
	{ 1.39730E+08, 3.49140E+06, 2.06579E+05, 5.65950E+04},
	{-2.56692E+08, 1.76440E+07, 6.26529E+05, 6.57440E+04},
	{ 0.00000E+00, 0.00000E+00, 9.30497E+05, 7.48920E+04}};


	init=0;
	ams=0.21580e1;
	anyu=795774.7155;

	bb=sqrt(bb);
	if(init==1)
	{
		goto bhodf211;
	}
	init=1;
bhodf211:;
	for(i=1;i<54;i++)
	{
		if(bb<bw[i])
		{
			goto bhodf212;
		}
	}
	*bh=anyu*(1.0-ams/bb);
	bb=bb*bb;
	*df=0.5/bb*(anyu-(*bh));
	return 0;
bhodf212:;
	i=i-1;
	b=bb-bw[i];
	*bh=((ch[i][0]*b+ch[i][1])*b+ch[i][2])*b+ch[i][3];
	*df=(3.0*ch[i][0]*b+2.0*ch[i][1])*b+ch[i][2];
	*bh=(*bh)/bb;
	bb=bb*bb;
	*df=0.5/bb*((*df)-(*bh));

	return 0;
}

int isus15c(double bb,double *bh,double *df)
{
	int init,i;
	double ams,anyu,b;
	double bw[38]={0.0000,
                   0.0931,
                   0.1950,
                   0.2929,
                   0.4101,
                   0.5085,
                   0.6072,
                   0.7160,
                   0.8120,
                   0.8978,
                   1.0001,
                   1.1013,
                   1.2055,
                   1.3034,
                   1.3932,      
                   1.5079,
                   1.6010,
                   1.6975,
                   1.8261,
                   1.9071,
                   2.0040,
                   2.0992,
                   2.1947,
                   2.2523,
                   2.2761,
                   2.2863,
                   2.3108,
                   2.3309,
                   2.3499,
                   2.3678,
                   2.3848,
                   2.4008,
                   2.4157,
                   2.4297,
                   2.4426,
                   2.4545,
                   2.4654,
				   2.4758};
	double ch[38][4]={{-5.04045E+04, 4.69363E+03, 1.10646E+03, 0.00000E+00},
	{ 1.70999E+04,-4.35081E+03, 6.69390E+02, 1.03032E+02},
	{-2.37210E+03, 1.68554E+02, 3.15300E+02, 1.44146E+02},
	{ 3.28752E+03,-5.61146E+02, 2.80109E+02, 1.74401E+02},
	{ 3.08881E+01, 4.31434E+02, 2.84040E+02, 2.04813E+02},
	{-1.84485E+03, 6.20938E+02, 3.69856E+02, 2.36975E+02},
	{-8.96265E+01, 3.98666E+02, 4.38510E+02, 2.77744E+02},
	{-4.61057E+03, 1.51584E+03, 5.22142E+02, 3.30103E+02},
	{-2.39085E+03, 1.46825E+03, 6.85701E+02, 3.90105E+02},
	{-5.39953E+03, 2.18981E+03, 8.84838E+02, 4.58230E+02},
	{-9.14276E+03, 2.90925E+03, 1.16335E+03, 5.65889E+02},
	{-6.65322E+03, 4.14425E+03, 1.47127E+03, 7.03854E+02},
	{-4.46476E+04, 1.21246E+04, 2.11859E+03, 8.94825E+02},
	{-5.29798E+04, 1.48006E+04, 3.20889E+03, 1.17646E+03},
	{-8.22662E+04, 3.01274E+04, 4.58536E+03, 1.54548E+03},
	{-1.76181E+05, 6.87349E+04, 8.24972E+03, 2.34372E+03},
	{-3.50993E+05, 1.10661E+05, 1.64683E+04, 3.56594E+03},
	{ 2.79031E+05, 3.79885E+04, 2.80189E+04, 5.86798E+03},
	{-1.25565E+06, 2.31840E+05, 5.16474E+04, 1.06954E+04},
	{-5.27020E+05, 1.32406E+05, 6.44916E+04, 1.57321E+04},
	{-2.03490E+06, 1.04959E+06, 7.53048E+04, 2.27480E+04},
	{ 5.70203E+06, 1.76201E+06, 2.19799E+05, 3.76688E+04},
	{-3.91107E+06, 6.18691E+05, 7.12058E+05, 7.96649E+04},
	{-2.56760E+06, 6.44291E+05, 7.44401E+05, 1.21997E+05},
	{-1.14603E+10, 1.56584E+08, 7.70752E+05, 1.40083E+05},
	{ 4.21654E+07,-1.70991E+06, 4.19388E+05, 1.52041E+05},
	{ 1.43382E+07,-2.74366E+05, 4.11763E+05, 1.61944E+05},
	{-1.27011E+07, 1.14635E+06, 4.18066E+05, 1.70210E+05},
	{-6.24747E+06, 8.70508E+05, 4.47872E+05, 1.78480E+05},
	{-4.84853E+06, 8.73037E+05, 4.73031E+05, 1.86740E+05},
	{-2.40343E+06, 1.18614E+06, 4.98511E+05, 1.95010E+05},
	{-2.45922E+07, 1.73643E+06, 5.34622E+05, 2.03280E+05},
	{ 4.46751E+05, 1.42336E+06, 5.69988E+05, 2.11550E+05},
	{-4.09535E+07, 2.92944E+06, 6.10105E+05, 2.19810E+05},
	{ 5.20625E+07, 1.87809E+06, 6.65239E+05, 2.28080E+05},
	{-2.68826E+07, 2.65463E+06, 7.32055E+05, 2.36350E+05},
	{-1.37235E+08, 2.85450E+06, 7.80344E+05, 2.44610E+05},
	{ 0.00000E+00, 0.00000E+00, 7.95188E+05, 2.52880E+05}};

	
	init=0;
	ams=0.21580e1;
	anyu=795774.7155;

	bb=sqrt(bb);
	if(init==1)
	{
		goto isus15c1;
	}
	init=1;
isus15c1:;
	for(i=1;i<54;i++)
	{
		if(bb<bw[i])
		{
			goto isus15c2;
		}
	}
	*bh=anyu*(1.0-ams/bb);
	bb=bb*bb;
	*df=0.5/bb*(anyu-(*bh));
	return 0;
isus15c2:;
	i=i-1;
	b=bb-bw[i];
	*bh=((ch[i][0]*b+ch[i][1])*b+ch[i][2])*b+ch[i][3];
	*df=(3.0*ch[i][0]*b+2.0*ch[i][1])*b+ch[i][2];
	*bh=(*bh)/bb;
	bb=bb*bb;
	*df=0.5/bb*((*df)-(*bh));

	return 0;
}

int src(double bb,double *bh,double *df)
{
	int init,i;
	double ams,anyu,b;
	double bw[34]={0.0000,
                   0.8000,
                   0.8127,
                   0.8276,
                   0.8445,
                   0.8628,
                   0.8816,
                   0.8997,
                   0.9163,
                   0.9308,
                   0.9432,
                   0.9535,
                   0.9620,
                   0.9690,
                   0.9747,    
                   0.9794,
                   0.9833,
                   0.9866,
                   0.9894,
                   0.9917,
                   0.9937,
                   0.9956,
                   0.9972,
                   0.9988,
                   1.0005,
                   1.0021,
                   1.0036,
                   1.0050,
                   1.0062,
                   1.0074,     
                   1.0084,
                   1.0094,
                   1.0102,
	               1.0109};

	double ch[34][4]={{-1.30316E+01, 1.34954E+01, 7.93319E+02, 0.00000E+00},
	{-4.61201E+04, 1.24551E+03, 7.89891E+02, 6.36620E+02},
	{-4.45527E+04, 1.86524E+03, 7.99211E+02, 6.46758E+02},
	{-6.52300E+04, 3.43893E+03, 8.25121E+02, 6.58933E+02},
	{-9.00284E+04, 5.58022E+03, 8.85464E+02, 6.73536E+02},
	{-1.35035E+05, 9.03939E+03, 9.99263E+02, 6.91067E+02},
	{-2.21334E+05, 1.48488E+04, 1.19582E+03, 7.12091E+02},
	{-3.92820E+05, 2.50785E+04, 1.51598E+03, 7.37333E+02},
	{-7.75391E+05, 4.39309E+04, 2.02385E+03, 7.67612E+02},
	{-1.54322E+06, 7.69395E+04, 2.80979E+03, 8.03971E+02},
	{-3.34577E+06, 1.38657E+05, 4.00485E+03, 8.47580E+02},
	{-6.91739E+06, 2.44635E+05, 5.79703E+03, 8.99942E+02},
	{-1.51605E+07, 4.40920E+05, 8.45784E+03, 9.62728E+02},
	{-3.16664E+07, 7.71206E+05, 1.23972E+04, 1.03809E+03},
	{-6.34479E+07, 1.33333E+06, 1.81161E+04, 1.12849E+03},
	{-1.35117E+08, 2.30130E+06, 2.64622E+04, 1.23703E+03},
	{-2.38648E+08, 3.70016E+06, 3.82469E+04, 1.36722E+03},
	{-4.80463E+08, 6.18972E+06, 5.47905E+04, 1.52351E+03},
	{-6.65679E+08, 9.15022E+06, 7.79334E+04, 1.71100E+03},
	{-1.37857E+09, 1.41788E+07, 1.09911E+05, 1.93604E+03},
	{-9.69806E+08, 1.69338E+07, 1.50434E+05, 2.20605E+03},
	{ 9.11949E+09, 6.16537E+06, 2.02436E+05, 2.53009E+03},
	{-5.87427E+08, 3.98813E+06, 2.96244E+05, 2.91890E+03},
	{ 1.09004E+09,-1.09177E+06, 3.04373E+05, 3.38546E+03},
	{-4.12441E+09, 2.03223E+07, 3.10474E+05, 3.91750E+03},
	{-1.88952E+09, 1.00568E+07, 3.43830E+05, 4.44940E+03},
	{ 2.26947E+10,-1.83759E+07, 3.61247E+05, 4.98140E+03},
	{-8.90683E+07, 1.98470E+05, 4.43238E+05, 5.51340E+03},
	{ 2.52615E+10,-3.03372E+07, 4.43330E+05, 6.04540E+03},
	{-5.24120E+10, 1.04791E+08, 4.79660E+05, 6.57740E+03},
	{ 1.26801E+08,-1.64867E+05, 5.32013E+05, 7.10940E+03},
	{-1.35072E+11, 2.74213E+08, 5.32064E+05, 7.64140E+03},
	{-7.55491E+09, 7.46306E+07, 7.11465E+05, 8.17340E+03},
	{ 0.00000E+00, 0.00000E+00, 8.04841E+05, 8.70540E+03}};

	init=0;
	ams=0.21580e1;
	anyu=795774.7155;

	bb=sqrt(bb);
	if(init==1)
	{
		goto src1;
	}
	init=1;
src1:;
	for(i=1;i<54;i++)
	{
		if(bb<bw[i])
		{
			goto src2;
		}
	}
	*bh=anyu*(1.0-ams/bb);
	bb=bb*bb;
	*df=0.5/bb*(anyu-(*bh));
	return 0;
src2:;
	i=i-1;
	b=bb-bw[i];
	*bh=((ch[i][0]*b+ch[i][1])*b+ch[i][2])*b+ch[i][3];
	*df=(3.0*ch[i][0]*b+2.0*ch[i][1])*b+ch[i][2];
	*bh=(*bh)/bb;
	bb=bb*bb;
	*df=0.5/bb*((*df)-(*bh));

	return 0;
}

int steel(double bb,double *bh,double *df)
{
	int init,i;
	double ams,anyu,b;

	double bw[32]={0.0000,
                   0.0490,
                   0.1010,
                   0.1500,
                   0.2000,
                   0.2990,
                   0.3990,
                   0.4990,
                   0.6010,
                   0.7000,
                   0.8010,
                   0.8990,
                   1.0010,
                   1.0990,
                   1.2000,     
                   1.3000,
                   1.4010,
                   1.4490,
                   1.5000,
                   1.5500,
                   1.6000,
                   1.6390,
                   1.6700,
                   1.7010,
                   1.7290,
                   1.7600,
                   1.7810,
                   1.8000,
                   1.8300,
                   1.8500,     
                   1.8750,
                   1.9000};

	double ch[32][4]={{-4.92643E+04,-1.01605E+04, 2.96309E+03,  .00000E+00},
	{ 1.06522E+04,-1.08536E+04, 1.61251E+03, 1.15000E+02},
	{-3.28171E+04, 3.84837E+02, 5.70141E+02, 1.71000E+02},
	{-1.13135E+04, 3.36207E+02, 3.71473E+02, 1.96000E+02},
	{ 2.39713E+03,-3.09152E+02, 3.20243E+02, 2.14000E+02},
	{-2.09717E+01, 1.06959E+02, 3.29514E+02, 2.45000E+02},
	{ 1.44659E+02, 1.82769E+02, 3.50276E+02, 2.79000E+02},
	{-1.41153E+03, 4.42004E+02, 3.91170E+02, 3.16000E+02},
	{-5.01029E+02, 3.26014E+02, 4.37282E+02, 3.59000E+02},
	{-2.37175E+03, 9.06425E+02, 4.87101E+02, 4.05000E+02},
	{-3.60831E+03, 1.23176E+03, 5.97616E+02, 4.61000E+02},
	{-3.62606E+03, 1.62150E+03, 7.35077E+02, 5.28000E+02},
	{-1.16315E+04, 3.49691E+03, 9.52686E+02, 6.16000E+02},
	{-1.59948E+04, 4.98784E+03, 1.30295E+03, 7.32000E+02},
	{-2.91278E+04, 1.03027E+04, 1.82101E+03, 8.98000E+02},
	{-3.20456E+04, 1.77666E+04, 3.00772E+03, 1.15400E+03},
	{-3.04884E+05, 5.34530E+04, 5.61588E+03, 1.60600E+03},
	{-1.56561E+05, 4.65695E+04, 8.64001E+03, 1.96500E+03},
	{-2.38904E+05, 8.25764E+04, 1.21685E+04, 2.50600E+03},
	{-3.57624E+04, 8.47011E+04, 1.86343E+04, 3.29100E+03},
	{ 2.18359E+05, 7.19486E+04, 2.68362E+04, 4.43000E+03},
	{-2.89934E+05, 7.37352E+04, 3.34446E+04, 5.59900E+03},
	{ 2.56506E+06,-1.04498E+03, 3.71803E+04, 6.69800E+03},
	{-1.65892E+06, 1.46839E+05, 4.45105E+04, 7.92600E+03},
	{ 1.89760E+05, 2.24365E+04, 4.88317E+04, 9.25100E+03},
	{ 8.08204E+06,-6.84244E+03, 5.07699E+04, 1.07920E+04},
	{-2.73398E+06, 8.98336E+04, 6.11751E+04, 1.19300E+04},
	{ 1.73982E+06,-5.86855E+04, 6.16279E+04, 1.31060E+04},
	{-9.25187E+06, 3.97327E+05, 6.28043E+04, 1.49490E+04},
	{ 3.76796E+06,-5.88007E+04, 6.75951E+04, 1.62900E+04},
	{-2.10923E+06, 2.87934E+05, 7.17200E+04, 1.80020E+04},
	{ 0.00000E+00, 0.00000E+00, 8.21619E+04, 1.99420E+04}};


	init=0;
	ams=0.21580e1;
	anyu=795774.7155;

	bb=sqrt(bb);
	if(init==1)
	{
		goto steel1;
	}
	init=1;
steel1:;
	for(i=1;i<54;i++)
	{
		if(bb<bw[i])
		{
			goto steel2;
		}
	}
	*bh=anyu*(1.0-ams/bb);
	bb=bb*bb;
	*df=0.5/bb*(anyu-(*bh));
	return 0;
steel2:;
	i=i-1;
	b=bb-bw[i];
	*bh=((ch[i][0]*b+ch[i][1])*b+ch[i][2])*b+ch[i][3];
	*df=(3.0*ch[i][0]*b+2.0*ch[i][1])*b+ch[i][2];
	*bh=(*bh)/bb;
	bb=bb*bb;
	*df=0.5/bb*((*df)-(*bh));

	return 0;
}

int ss400(double bb,double *bh,double *df)
{
	int init,i;
	double ams,anyu,b;

	double bw[50]={0.00000000,
               0.04900000,
               0.10100000,
               0.15000001,
               0.20000000,
               0.29899999,
               0.39899999,
               0.49900001,
               0.60100001,
               0.69999999,
               0.80100000,
               0.89899999,
               1.00100005,
               1.09899998,
               1.20000005,      
               1.29999995,
               1.40100002,
               1.44900000,
               1.50000000,
               1.54999995,
               1.60000002,
               1.63900006,
               1.66999996,
               1.70099998,
               1.72899997,
               1.75999999,
               1.78100002,
               1.79999995,
               1.83000004,
               1.85000002,      
               1.87500000,
               1.89999998,
               1.93589997,
               1.96990001,
               2.00189996,
               2.03189993,
               2.05999994,
               2.08610010,
               2.11030006,
               2.13249993,
               2.15280008,
               2.17100000,
               2.18740010,
               2.20180011,
               2.21420002,     
               2.22460008,
               2.23309994,
               2.23970008,
               2.24419999,
               2.24780011};

	double ch[50][4]={{-3.05885E+05, 1.49884E+04, 2.34694E+03, 0.00000E+00},
	{ 1.06522E+04,-1.08536E+04, 1.61251E+03, 1.15000E+02},
	{-3.28171E+04, 3.84837E+02, 5.70141E+02, 1.71000E+02},
	{-1.13135E+04, 3.36207E+02, 3.71473E+02, 1.96000E+02},
	{ 2.39713E+03,-3.09152E+02, 3.20243E+02, 2.14000E+02},
	{-2.09717E+01, 1.06959E+02, 3.29514E+02, 2.45000E+02},
	{ 1.44659E+02, 1.82769E+02, 3.50276E+02, 2.79000E+02},
	{-1.41153E+03, 4.42004E+02, 3.91170E+02, 3.16000E+02},
	{-5.01029E+02, 3.26014E+02, 4.37282E+02, 3.59000E+02},
	{-2.37175E+03, 9.06425E+02, 4.87101E+02, 4.05000E+02},
	{-3.60831E+03, 1.23176E+03, 5.97616E+02, 4.61000E+02},
	{-3.62606E+03, 1.62150E+03, 7.35077E+02, 5.28000E+02},
	{-1.16315E+04, 3.49691E+03, 9.52686E+02, 6.16000E+02},
	{-1.59948E+04, 4.98784E+03, 1.30295E+03, 7.32000E+02},
	{-2.91278E+04, 1.03027E+04, 1.82101E+03, 8.98000E+02},
	{-3.20456E+04, 1.77666E+04, 3.00772E+03, 1.15400E+03},
	{-3.04884E+05, 5.34530E+04, 5.61588E+03, 1.60600E+03},
	{-1.56561E+05, 4.65695E+04, 8.64001E+03, 1.96500E+03},
	{-2.38904E+05, 8.25764E+04, 1.21685E+04, 2.50600E+03},
	{-3.57624E+04, 8.47011E+04, 1.86343E+04, 3.29100E+03},
	{ 2.18359E+05, 7.19486E+04, 2.68362E+04, 4.43000E+03},
	{-2.89934E+05, 7.37352E+04, 3.34446E+04, 5.59900E+03},
	{ 2.56506E+06,-1.04498E+03, 3.71803E+04, 6.69800E+03},
	{-1.65892E+06, 1.46839E+05, 4.45105E+04, 7.92600E+03},
	{ 1.89760E+05, 2.24365E+04, 4.88317E+04, 9.25100E+03},
	{ 8.08204E+06,-6.84244E+03, 5.07699E+04, 1.07920E+04},
	{-2.73398E+06, 8.98336E+04, 6.11751E+04, 1.19300E+04},
	{ 1.73982E+06,-5.86855E+04, 6.16279E+04, 1.31060E+04},
	{-9.25187E+06, 3.97327E+05, 6.28043E+04, 1.49490E+04},
	{ 4.41396E+06,-7.49504E+04, 6.75951E+04, 1.62900E+04},
	{-6.45714E+06, 3.80482E+05, 7.21237E+04, 1.80020E+04},
	{ 4.63436E+05, 3.10128E+03, 7.90407E+04, 1.99420E+04},
	{-8.13620E+05, 1.21193E+05, 8.10552E+04, 2.28050E+04},
	{-2.93412E+05, 1.02959E+05, 8.64746E+04, 2.56690E+04},
	{-5.72171E+05, 1.27301E+05, 9.21627E+04, 2.85320E+04},
	{-3.87078E+05, 1.40067E+05, 9.82558E+04, 3.13960E+04},
	{-1.19650E+06, 2.02959E+05, 1.05211E+05, 3.42590E+04},
	{-9.62705E+05, 2.29388E+05, 1.13360E+05, 3.71220E+04},
	{-2.84223E+06, 3.42102E+05, 1.22771E+05, 3.99860E+04},
	{-1.53567E+06, 3.92004E+05, 1.33758E+05, 4.28490E+04},
	{-7.92314E+06, 6.68032E+05, 1.47775E+05, 4.57130E+04},
	{-6.81541E+06, 7.43143E+05, 1.64218E+05, 4.85760E+04},
	{-2.16387E+07, 1.40849E+06, 1.83094E+05, 5.14390E+04},
	{-3.54771E+07, 2.10861E+06, 2.10197E+05, 5.43030E+04},
	{-9.17924E+07, 3.76786E+06, 2.46126E+05, 5.71660E+04},
	{-2.70768E+08, 7.25648E+06, 2.94712E+05, 6.00300E+04},
	{ 5.59635E+07, 1.09024E+07, 3.59383E+05, 6.28930E+04},
	{-1.80453E+09, 3.60864E+07, 5.10612E+05, 6.57560E+04},
	{-5.36146E+09, 3.86038E+07, 7.25762E+05, 6.86200E+04},
	{ 0.00000E+00, 0.00000E+00, 7.95775E+05, 7.14830E+04}};


	init=0;
	ams=0.21580e1;
	anyu=795774.7155;

	bb=sqrt(bb);
	if(init==1)
	{
		goto ss4001;
	}
	init=1;
ss4001:;
	for(i=1;i<54;i++)
	{
		if(bb<bw[i])
		{
			goto ss4002;
		}
	}
	*bh=anyu*(1.0-ams/bb);
	bb=bb*bb;
	*df=0.5/bb*(anyu-(*bh));
	return 0;
ss4002:;
	i=i-1;
	b=bb-bw[i];
	*bh=((ch[i][0]*b+ch[i][1])*b+ch[i][2])*b+ch[i][3];
	*df=(3.0*ch[i][0]*b+2.0*ch[i][1])*b+ch[i][2];
	*bh=(*bh)/bb;
	bb=bb*bb;
	*df=0.5/bb*((*df)-(*bh));

	return 0;
}







