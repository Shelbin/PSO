//����ȺPSO�㷨

#include<stdio.h>
#include<math.h>
#include<time.h>
#include <sys/time.h> 
#include<stdlib.h>
#include<mpi.h>

#define	PI 3.141592653589	
#define P_num 20000          			//������Ŀ
#define dim 50					       //ά��
#define N   10						  // ÿ������ִ����Ŀ

#define low -100           			  //������Χ
#define high 100

#define iter_num 1000
#define V_max 20          		//�ٶȷ�Χ

#define c1 2					//�����ٶ�ʱ�ļ��ٳ���
#define c2 2					//�����ٶ�ʱ�ļ��ٳ���
#define w 0.5					//�����ٶ�ʱ��Ȩ��
#define alp 1					//����λ��ʱ��Ȩ��

double particle[P_num][dim];           //����λ�ü���
double particle_loc_best[P_num][dim];  //ÿ������ֲ���������
double particle_loc_fit[P_num];        //����ľֲ�������Ӧ��,�оֲ����������������
double particle_glo_best[dim];         //ȫ����������
double gfit;                           //ȫ��������Ӧ��,��ȫ�����������������
double particle_v[P_num][dim];         //��¼ÿ������ĵ�ǰ���ٶ�����
double particle_fit[P_num];            //��¼ÿ�����ӵĵ�ǰ����Ӧ��

/*ʱ�亯��*/
double mytime()
{
    double ts = 0.0;
    struct timeval mt;
    gettimeofday(&mt,(struct timezone*)0);
    ts = (double)(mt.tv_sec+mt.tv_usec*1.0e-6);
    return (ts);
}

int get_rand(int comm_sz){
  srand((unsigned)time(NULL));
  return rand()%(P_num/comm_sz);
}

double Sphere(double a[])
{

	int i;
	double sum=0.0;
	for(i=0; i<dim; i++)
	{
		sum+=a[i]*a[i];
	}
	return sum;
}
double Rosenbrock(double a[])
{

	int i;
	double sum=0.0;
	for(i=0;i<dim-1; i++)
	{
        sum+= 100*(a[i+1]-a[i]*a[i])*(a[i+1]-a[i]*a[i])+(a[i]-1)*(a[i]-1);
	}
  	    return sum;
}

double Rastrigin(double a[])
{
	int i;
	double sum=0.0;
	for(i=0;i<dim;i++)
	{
		sum+=a[i]*a[i]-10.0*cos(2*PI*a[i])+10.0;
	}
	return sum;
}

double fitness(double a[])             //��Ӧ�Ⱥ���
{
	return Rastrigin(a);
}

/*��ʼ�����ӣ�������ʼ���ٶȣ�λ�ã����Ե���Ӧ���Լ�*/
void initial()
{
	int i,j;
	for(i=0; i<P_num; i++)            			//�����������
	{
		for(j=0; j<dim; j++)
		{
			particle[i][j] = low+(high-low)*1.0*rand()/RAND_MAX;    //��ʼ��Ⱥ��,λ�����
			particle_loc_best[i][j] = particle[i][j];               //����ǰ���Ž��д��ֲ����ż���
			particle_v[i][j] = -V_max+2*V_max*1.0*rand()/RAND_MAX;    //��ʼ���ٶ�
		}
	}
	
	for(i=0; i<P_num; i++)            			//����ÿ�����ӵ���Ӧ��
	{
		particle_fit[i] = fitness(particle[i]);
		particle_loc_fit[i] = particle_fit[i];			//��ǰ�ֲ���Ӧ��
	}
	
	gfit = particle_loc_fit[0];      		//��ʼ��ȫ��������Ӧ��
	
	j=0;
	
	for(i=1; i<P_num; i++)
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;							//��¼�±�
			
		}
	}
	
	for(i=0; i<dim; i++)             		//����ȫ����������  ά�ȿ������Ϊ��Ⱥ����
	{
		particle_glo_best[i] = particle_loc_best[j][i];
	}
}

/*����λ�ú���*/
void renew_particle(int a,int b)
{
	int i,j;
	for(i=a; i<b; i++)            //���¸���λ������λ��
	{
		for(j=0; j<dim; j++)
		{
			particle[i][j] +=  alp*particle_v[i][j];    
			if(particle[i][j] > high)
			{
				particle[i][j] = high;
			}
			if(particle[i][j] < low)
			{
				particle[i][j] = low;
			}
		}
	}
}

/*���������ٶ�*/
void renew_var(int a,int b,int *m)
{
	int i, j;
	for(i=a;i<b;i++)            						//����ÿ�����ӵ���Ӧ��
	{
		particle_fit[i] = fitness(particle[i]);			//�����������Ӧ��
		
		if(particle_fit[i] < particle_loc_fit[i])      //���¸���ֲ�����ֵ
		{
			particle_loc_fit[i] = particle_fit[i];
			
			for(j=0; j<dim; j++)      				 // ���¾ֲ���������
			{
				particle_loc_best[i][j] = particle[i][j];
			}
		}
	}
	
	for(i=a,j=-1; i<b; i++)                   	//����ȫ�ֱ���
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;
		
		}
	}
	
	
	if(j != -1)
	{
		for(i=0; i<dim; i++)             //����ȫ����������
		{
			particle_glo_best[i] = particle_loc_best[j][i];
		}
	}
	*m = j;
	
	for(i=a; i<b; i++)    //���¸����ٶ�
	{
		for(j=0; j<dim; j++)
		{
			particle_v[i][j]=w*particle_v[i][j]+
				c1*1.0*rand()/RAND_MAX*(particle_loc_best[i][j]-particle[i][j])+
				c2*1.0*rand()/RAND_MAX*(particle_glo_best[j]-particle[i][j]);
			if(particle_v[i][j] > V_max)
			{
				particle_v[i][j] = V_max;
			}
			if(particle_v[i][j] < -V_max)
			{
				particle_v[i][j] = -V_max;
			}
		}
	}
}


int main()
{
 
  
  double tstart,tstop;
	//freopen("result.txt","a+",stdout);
  int i=0,j;
  int a,b;
  double *buf3;
  int *buf4;
  int n,k,m,ii;
  
  
  int my_rank,comm_sz;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
  
  //if(my_rank==0){
     tstart = mytime(); 
  	 srand((unsigned)time(NULL));			//������������ĳ�ʼ������
  	 initial();						      //��ʼ��Ⱥ��
  /* MPI_Bcast(particle,1000000,MPI_DOUBLE,0,MPI_COMM_WORLD); 
     MPI_Bcast(particle_loc_best,1000000,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_loc_fit,P_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_glo_best,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);               
     MPI_Bcast(&gfit,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_v,1000000,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_fit,P_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }*/
  
 
	   
 
   //MPI_Barrier(MPI_COMM_WORLD);          //ͬ��
   
   if(my_rank==0){
     buf3 = (double *)malloc(comm_sz*sizeof(double));
     buf4 = (int *)malloc(comm_sz*sizeof(int));
   }
 
   MPI_Barrier(MPI_COMM_WORLD);
   
   a=my_rank*(P_num/comm_sz);
   b=a+(P_num/comm_sz);
   
  
  do{
	 
    for(j=0;j<N;j++){
        renew_particle(a,b);
        renew_var(a,b,&m);
    }
    
    MPI_Gather(&m,1,MPI_INT,buf4,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&gfit,1,MPI_DOUBLE,buf3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
    
   
    if(my_rank==0){
    
      gfit=buf3[0]; 
      for(i=1,j=0; i<comm_sz; i++)                   //����ȫ�ֱ���
      {
        if(buf3[i]<gfit){
         gfit = buf3[i];
          j = i;
        }
      }
      
     printf("gfit=%lf\n",gfit);
      for(i=0; i<dim; i++)             //����ȫ����������
      {
        particle_glo_best[i] = particle_loc_best[buf4[j]][i];
      }
    
    }
    
	
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=0;i<50;i++){
	   particle[get_rand(comm_sz)+a][i]=particle_glo_best[i];
	}
	
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank==0){
    
	    MPI_Bcast(particle_glo_best,50,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&gfit,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	
    MPI_Barrier(MPI_COMM_WORLD);
  }while(gfit>20);
  
	tstop=mytime();
	printf("The number of granule:%d\n",P_num);
	printf("Dimensionality:%d\n",dim);
	printf("Best result:%.10lf\n", gfit);
	printf("Spend time=%lfs\n",tstop-tstart);
	
	return 0;
}
