//����ȺPSO�㷨

#include<stdio.h>
#include<math.h>
#include<time.h>
#include <sys/time.h> 
#include<stdlib.h>

#define	PI 3.141592653589	
#define P_num 20000          		//������Ŀ
#define dim 50					//ά��

#define low -100           		//������Χ
#define high 100
#define V_max 20          		//�ٶȷ�Χ

#define c1 2
#define c2 2
#define w 0.5
#define alp 1
#define N   10	

double particle[P_num][dim];           //����λ�ü���
double particle_loc_best[P_num][dim];  //ÿ������ֲ���������
double particle_loc_fit[P_num];        //����ľֲ�������Ӧ��,�оֲ����������������
double particle_glo_best[dim];         //ȫ����������
double gfit;                           //ȫ��������Ӧ��,��ȫ�����������������
double particle_v[P_num][dim];         //��¼ÿ������ĵ�ǰ���ٶ�����
double particle_fit[P_num];            //��¼ÿ�����ӵĵ�ǰ����Ӧ��


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
double fitness(double a[])             //��Ӧ�Ⱥ���
{
	return Rastrigin(a);
}
void initial()
{
	int i,j;
	for(i=0; i<P_num; i++)            	//�����������
	{
		for(j=0; j<dim; j++)
		{
			particle[i][j] = low+(high-low)*1.0*rand()/RAND_MAX;    //��ʼ��Ⱥ��
			particle_loc_best[i][j] = particle[i][j];               //����ǰ���Ž��д��ֲ����ż���
			particle_v[i][j] = -V_max+2*V_max*1.0*rand()/RAND_MAX;    //�ٶ�
		}
	}
	
	for(i=0; i<P_num; i++)            //����ÿ�����ӵ���Ӧ��
	{
		particle_fit[i] = fitness(particle[i]);
		particle_loc_fit[i] = particle_fit[i];
	}
	
	gfit = particle_loc_fit[0];      //�ҳ�ȫ������
	j=0;
	
	for(i=1; i<P_num; i++)
	{
		if(particle_loc_fit[i]<gfit)
		{
			gfit = particle_loc_fit[i];
			j = i;
		}
	}
	
	for(i=0; i<dim; i++)             //����ȫ����������
	{
		particle_glo_best[i] = particle_loc_best[j][i];
	}
}

void renew_particle()
{
	int i,j;
	for(i=0; i<P_num; i++)            //���¸���λ������λ��
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
	for(i=0; i<P_num; i++)            					//����ÿ�����ӵ���Ӧ��
	{
		particle_fit[i] = fitness(particle[i]);
		if(particle_fit[i] < particle_loc_fit[i])      //���¸���ֲ�����ֵ
		{
			particle_loc_fit[i] = particle_fit[i];
			
			for(j=0; j<dim; j++)       // ���¾ֲ���������
			{
				particle_loc_best[i][j] = particle[i][j];
			}
		}
	}
	
	for(i=0,j=-1; i<P_num; i++)                   //����ȫ�ֱ���
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
		for(i=0; i<dim; i++)             //����ȫ����������
		{
			particle_glo_best[i] = particle_loc_best[j][i];
		}
	}
	
	for(i=0; i<P_num; i++)    //���¸����ٶ�
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
	initial();						//��ʼ��Ⱥ��
	
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
