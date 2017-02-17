//粒子群PSO算法

#include<stdio.h>
#include<math.h>
#include<time.h>
#include <sys/time.h> 
#include<stdlib.h>
#include<mpi.h>

#define	PI 3.141592653589	
#define P_num 20000          			//粒子数目
#define dim 50					       //维度
#define N   10						  // 每条进程执行数目

#define low -100           			  //搜索域范围
#define high 100

#define iter_num 1000
#define V_max 20          		//速度范围

#define c1 2					//更新速度时的加速常量
#define c2 2					//更新速度时的加速常量
#define w 0.5					//更新速度时的权重
#define alp 1					//更新位置时的权重

double particle[P_num][dim];           //个体位置集合
double particle_loc_best[P_num][dim];  //每个个体局部最优向量
double particle_loc_fit[P_num];        //个体的局部最优适应度,有局部最优向量计算而来
double particle_glo_best[dim];         //全局最优向量
double gfit;                           //全局最优适应度,有全局最优向量计算而来
double particle_v[P_num][dim];         //记录每个个体的当前代速度向量
double particle_fit[P_num];            //记录每个粒子的当前代适应度

/*时间函数*/
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

double fitness(double a[])             //适应度函数
{
	return Rastrigin(a);
}

/*初始化粒子，包括初始化速度，位置，各自的适应度以及*/
void initial()
{
	int i,j;
	for(i=0; i<P_num; i++)            			//随机生成粒子
	{
		for(j=0; j<dim; j++)
		{
			particle[i][j] = low+(high-low)*1.0*rand()/RAND_MAX;    //初始化群体,位置随机
			particle_loc_best[i][j] = particle[i][j];               //将当前最优结果写入局部最优集合
			particle_v[i][j] = -V_max+2*V_max*1.0*rand()/RAND_MAX;    //初始化速度
		}
	}
	
	for(i=0; i<P_num; i++)            			//计算每个粒子的适应度
	{
		particle_fit[i] = fitness(particle[i]);
		particle_loc_fit[i] = particle_fit[i];			//当前局部适应度
	}
	
	gfit = particle_loc_fit[0];      		//初始化全局最优适应度
	
	j=0;
	
	for(i=1; i<P_num; i++)
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;							//记录下标
			
		}
	}
	
	for(i=0; i<dim; i++)             		//更新全局最优向量  维度可以理解为种群个数
	{
		particle_glo_best[i] = particle_loc_best[j][i];
	}
}

/*更新位置函数*/
void renew_particle(int a,int b)
{
	int i,j;
	for(i=a; i<b; i++)            //更新个体位置生成位置
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

/*更新粒子速度*/
void renew_var(int a,int b,int *m)
{
	int i, j;
	for(i=a;i<b;i++)            						//计算每个粒子的适应度
	{
		particle_fit[i] = fitness(particle[i]);			//计算该粒子适应度
		
		if(particle_fit[i] < particle_loc_fit[i])      //更新个体局部最优值
		{
			particle_loc_fit[i] = particle_fit[i];
			
			for(j=0; j<dim; j++)      				 // 更新局部最优向量
			{
				particle_loc_best[i][j] = particle[i][j];
			}
		}
	}
	
	for(i=a,j=-1; i<b; i++)                   	//更新全局变量
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;
		
		}
	}
	
	
	if(j != -1)
	{
		for(i=0; i<dim; i++)             //更新全局最优向量
		{
			particle_glo_best[i] = particle_loc_best[j][i];
		}
	}
	*m = j;
	
	for(i=a; i<b; i++)    //更新个体速度
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
  	 srand((unsigned)time(NULL));			//随机数发生器的初始化函数
  	 initial();						      //初始化群体
  /* MPI_Bcast(particle,1000000,MPI_DOUBLE,0,MPI_COMM_WORLD); 
     MPI_Bcast(particle_loc_best,1000000,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_loc_fit,P_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_glo_best,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);               
     MPI_Bcast(&gfit,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_v,1000000,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(particle_fit,P_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }*/
  
 
	   
 
   //MPI_Barrier(MPI_COMM_WORLD);          //同步
   
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
      for(i=1,j=0; i<comm_sz; i++)                   //更新全局变量
      {
        if(buf3[i]<gfit){
         gfit = buf3[i];
          j = i;
        }
      }
      
     printf("gfit=%lf\n",gfit);
      for(i=0; i<dim; i++)             //更新全局最优向量
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
