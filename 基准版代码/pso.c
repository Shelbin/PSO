//粒子群PSO算法

#include<stdio.h>
#include<math.h>
#include<time.h>
#include <sys/time.h> 
#include<stdlib.h>

#define	PI 3.141592653589	
#define P_num 20000          		//粒子数目
#define dim 50					//维度

#define low -100           		//搜索域范围
#define high 100
#define V_max 20          		//速度范围

#define c1 2
#define c2 2
#define w 0.5
#define alp 1
#define N   10	

double particle[P_num][dim];           //个体位置集合
double particle_loc_best[P_num][dim];  //每个个体局部最优向量
double particle_loc_fit[P_num];        //个体的局部最优适应度,有局部最优向量计算而来
double particle_glo_best[dim];         //全局最优向量
double gfit;                           //全局最优适应度,有全局最优向量计算而来
double particle_v[P_num][dim];         //记录每个个体的当前代速度向量
double particle_fit[P_num];            //记录每个粒子的当前代适应度


double mytime()
{
    double ts = 0.0;
    struct timeval mt;
    gettimeofday(&mt,(struct timezone*)0);
    ts = (double)(mt.tv_sec+mt.tv_usec*1.0e-6);
    return (ts);
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
void initial()
{
	int i,j;
	for(i=0; i<P_num; i++)            	//随机生成粒子
	{
		for(j=0; j<dim; j++)
		{
			particle[i][j] = low+(high-low)*1.0*rand()/RAND_MAX;    //初始化群体
			particle_loc_best[i][j] = particle[i][j];               //将当前最优结果写入局部最优集合
			particle_v[i][j] = -V_max+2*V_max*1.0*rand()/RAND_MAX;    //速度
		}
	}
	
	for(i=0; i<P_num; i++)            //计算每个粒子的适应度
	{
		particle_fit[i] = fitness(particle[i]);
		particle_loc_fit[i] = particle_fit[i];
	}
	
	gfit = particle_loc_fit[0];      //找出全局最优
	j=0;
	
	for(i=1; i<P_num; i++)
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;
		}
	}
	
	for(i=0; i<dim; i++)             //更新全局最优向量
	{
		particle_glo_best[i] = particle_loc_best[j][i];
	}
}

void renew_particle()
{
	int i,j;
	for(i=0; i<P_num; i++)            //更新个体位置生成位置
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

void renew_var(int *m)
{
	int i, j;
	for(i=0; i<P_num; i++)            					//计算每个粒子的适应度
	{
		particle_fit[i] = fitness(particle[i]);
		if(particle_fit[i] < particle_loc_fit[i])      //更新个体局部最优值
		{
			particle_loc_fit[i] = particle_fit[i];
			
			for(j=0; j<dim; j++)       // 更新局部最优向量
			{
				particle_loc_best[i][j] = particle[i][j];
			}
		}
	}
	
	for(i=0,j=-1; i<P_num; i++)                   //更新全局变量
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;
		}
	}
	*m = gfit;
	printf("%lf\n",gfit);
	if(j != -1)
	{
		for(i=0; i<dim; i++)             //更新全局最优向量
		{
			particle_glo_best[i] = particle_loc_best[j][i];
		}
	}
	
	for(i=0; i<P_num; i++)    //更新个体速度
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
	int i=0;
	int gfit;
	tstart = mytime();
	srand((unsigned)time(NULL));
	initial();						//初始化群体
	
	do{
		for(i=0;i<N;i++){
		   renew_particle();
		   renew_var(&gfit);
		}
	}while(gfit>16);
	
	tstop=mytime();
	printf("The number of granule:%d\n",P_num);
	printf("Dimensionality:%d\n",dim);
	printf("Best result:%lf\n", gfit);
	printf("Spend time=%lfs\n",tstop-tstart);
	
	return 0;
}
