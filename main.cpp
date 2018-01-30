//
//  main.cpp
//  telomere-length shortest to shortest sim
//
//  Created by Ye yusong on 17/9/4.
//  Copyright © 2017年 Ye yusong. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <time.h>

//geometric distribution generater

double geg(double cumpro[3000])
{
    
    double rd=((double)rand()/(RAND_MAX));
    double ksi=0;
    
    for(int i=0;i<3000;i++)
    {if(rd<cumpro[i])
    {ksi=(i+1);
        break;}
    }
    if(rd>cumpro[2999])
        ksi=3000;
    
    
    return ksi;
}

double SampleNormal()
{
    double u = ((double)rand()/(RAND_MAX)) * 2 - 1;
    double v = ((double)rand()/(RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return SampleNormal();
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}


//main function

int main(int argc, const char * argv[])
{
    
    srand(time(0));
    rand();
    rand();
    
    //parameter
    
    double l0=0;
    double lmax=20000;
    double u=230;
    double M=1;
    double kr=0.1;
    double alpha=0.001;
    
    
    //simulation samples
    int samples=5000;
    //simulation loops
    int loops=10;
    
    //simulation pd
    int maxpd=500;
    int pd1=10;
    int pd2=33;
    int pd3=49;
    int pd4=55;
    int pd5=100;
    int pd6=200;
    int pd7=300;
    int pd8=400;
    int pd9=500;
    
    
    //range
    double std=2500;
    double mean=8000;
    
    
    // cells telomere
    
    double telomother[samples][46];
    
    double telodaughter1[samples][46];
    int labeldaughter1[samples];
    for(int i=0;i<samples;i++)
    {labeldaughter1[i]=0;}
    
    double telodaughter2[samples][46];
    int labeldaughter2[samples];
    for(int i=0;i<samples;i++)
    {labeldaughter2[i]=0;}
    
    double daughtermin1;
    double daughtermin2;
    double pol;
    
    
    
    
    //计算几何分布
    double p=0.005;
    double progeo[3000];
    double cumpro[3000];
    for(int i=0;i<3000;i++)
    {cumpro[i]=0;}
    for(int i=0;i<3000;i++)
    {progeo[i]=p*(pow((1-p),(double)(i+1)));}
    for(int i=0;i<3000;i++)
    {
        for(int j=0;j<i+1;j++)
        {
            cumpro[i]=cumpro[i]+progeo[j];
        }
    }
    
    //statisitical number
    double population[maxpd+1];
    
    double freq1[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq1[i]=0;}
    
    double freq2[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq2[i]=0;}
    
    double freq3[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq3[i]=0;}
    
    double freq4[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq4[i]=0;}
    
    double freq5[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq1[i]=0;}
    
    double freq6[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq1[i]=0;}
    
    double freq7[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq1[i]=0;}
    
    double freq8[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq1[i]=0;}
    
    double freq9[(int)(lmax+1)];
    for(int i=0;i<(lmax+1);i++)
    {freq1[i]=0;}
    
    
    using namespace std;
    ofstream file;
    file.open("result.txt");
    file<<endl;
    
    
    // simulation times
    
    
    for(int loop=0;loop<loops;loop++)
    {
        
    //initial mother cell
    for(int i=0;i<samples;i++)
        {
            for(int j=0;j<46;j++)
            {telomother[i][j]=double (floor((SampleNormal()*std+mean)));}
        }
        
        file<<endl<<(loop+1)<<endl;
        
    //initial population
        population[0]=1;
        for(int pop=1;pop<(maxpd+1);pop++)
        {population[pop]=0;}
    
    //time pd
    for(int i=1;i<(maxpd+1);i++)
    {
        //重新标记死亡标签(imp)
        for(int lab=0;lab<samples;lab++)
        {labeldaughter1[lab]=0;}
        
        for(int lab=0;lab<samples;lab++)
        {labeldaughter2[lab]=0;}
        
        
     // cell number
        for(int k=0;k<(samples);k++)
        {
            //redefine min
            daughtermin1=5000;
            daughtermin2=5000;
            
            //telomere number
            // first shortening
            
            for(int j=0;j<46;j++)
            {
                
                double rd0=((double)rand()/(RAND_MAX));
                
                if(rd0<0.5)
                {telodaughter1[k][j]=telomother[k][j]-u;
                    telodaughter2[k][j]=telomother[k][j];}
                else
                {telodaughter1[k][j]=telomother[k][j];
                    telodaughter2[k][j]=telomother[k][j]-u;}
            }
            
            
                
                //death judging
                
                for(int z=0;z<46;z++)
                {
                    if(daughtermin1>telodaughter1[k][z])
                    {daughtermin1=telodaughter1[k][z];}
                }
                for(int z=0;z<46;z++)
                {
                    if(daughtermin2>telodaughter2[k][z])
                    {daughtermin2=telodaughter2[k][z];}
                }
                
                if(daughtermin1<0)
                {labeldaughter1[k]=1;}
                else
                {double rdz=((double)rand()/(RAND_MAX));
                    if(rdz<(double (exp(-alpha*daughtermin1))))
                    {labeldaughter1[k]=1;}
                }

                if(daughtermin2<0)
                {labeldaughter2[k]=1;}
                else
                {double rdz=((double)rand()/(RAND_MAX));
                    if(rdz<(double (exp(-alpha*daughtermin2))))
                    {labeldaughter2[k]=1;}
                }
                
            //lengthing
                
                if(labeldaughter1[k]==0)
                {for(int j=0;j<46;j++)
                {
                    double pol=(2*M)/(1+2*M+pow((1+4*kr*telodaughter1[k][j]*(1+M)),0.5));
                    
                    
                    double rd1=((double)rand()/(RAND_MAX));
                    double ksi1=geg(cumpro);
                    
                    if(rd1<pol)
                    {telodaughter1[k][j]=telodaughter1[k][j]+ksi1;}
                    if(telodaughter1[k][j]>lmax)
                    {telodaughter1[k][j]=lmax;}
                }
                }
                if(labeldaughter2[k]==0)
                {for(int j=0;j<46;j++)
                {
                    double pol=(2*M)/(1+2*M+pow((1+4*kr*telodaughter2[k][j]*(1+M)),0.5));
                    
                    
                    double rd1=((double)rand()/(RAND_MAX));
                    double ksi1=geg(cumpro);
                    
                    if(rd1<pol)
                    {telodaughter2[k][j]=telodaughter2[k][j]+ksi1;}
                    if(telodaughter2[k][j]>lmax)
                    {telodaughter2[k][j]=lmax;}
                }
                }
            
                
            }
        
        
        
        
        
        //judging population deathing
        int label1=0;
        for(int j=0;j<samples;j++)
        {label1=label1+labeldaughter1[j];}
        if(label1==samples)
        {break;}
        
    
        
        int label2=0;
        for(int j=0;j<samples;j++)
        {label2=label2+labeldaughter2[j];}
        if(label1==samples)
        {break;}
        
        double d=(((double)label1+(double)label2)/10000);
        population[i]=population[i-1]*2*(1-d);
        
        file<<d<<"\t";
        
        if(population[i]<1)
        {break;}
        
            //initial mother cells
        
        for(int j=0;j<samples;j++)
        {
            int rd1=0;
            double rd2=(double)(rand()/(RAND_MAX));
            if(rd2<0.5)
            {
                
            do{
                rd1=((int)floor(((double)rand()/(RAND_MAX))*samples));
            {
            for(int k=0;k<46;k++)
            {telomother[j][k]=telodaughter1[rd1][k];}}}
                while(labeldaughter1[rd1]==1);
            }
            else
            {
                
                do{
                    rd1=((int)floor(((double)rand()/(RAND_MAX))*samples));
                    {
                        for(int k=0;k<46;k++)
                        {telomother[j][k]=telodaughter1[rd1][k];}}}
                while(labeldaughter1[rd1]==1);
            }
        }
        
            
            //statistical
            
            if(i==pd1)
            {for(int k=0;k<samples;k++)
            {
            {for(int j=0;j<46;j++)
            {freq1[(int)(telomother[k][j])]=freq1[(int)(telomother[k][j])]+1;}}}}
        
            if(i==pd2)
            {for(int k=0;k<samples;k++)
            {
                {for(int j=0;j<46;j++)
                {freq2[(int)(telomother[k][j])]=freq2[(int)(telomother[k][j])]+1;}}
            }}
            
            if(i==pd3)
            {for(int k=0;k<samples;k++)
            {
                {for(int j=0;j<46;j++)
                {freq3[(int)(telomother[k][j])]=freq3[(int)(telomother[k][j])]+1;}}
            }}
        
            if(i==pd4)
            {for(int k=0;k<samples;k++)
            {
                {for(int j=0;j<46;j++)
                {freq4[(int)(telomother[k][j])]=freq4[(int)(telomother[k][j])]+1;}}
            }}
        
        if(i==pd5)
        {for(int k=0;k<samples;k++)
        {
            {for(int j=0;j<46;j++)
            {freq5[(int)(telomother[k][j])]=freq5[(int)(telomother[k][j])]+1;}}
        }}
        
        if(i==pd6)
        {for(int k=0;k<samples;k++)
        {
            {for(int j=0;j<46;j++)
            {freq6[(int)(telomother[k][j])]=freq6[(int)(telomother[k][j])]+1;}}
        }}
        
        if(i==pd7)
        {for(int k=0;k<samples;k++)
        {
            {for(int j=0;j<46;j++)
            {freq7[(int)(telomother[k][j])]=freq7[(int)(telomother[k][j])]+1;}}
        }}
        
        if(i==pd8)
        {for(int k=0;k<samples;k++)
        {
            {for(int j=0;j<46;j++)
            {freq8[(int)(telomother[k][j])]=freq8[(int)(telomother[k][j])]+1;}}
        }}
        
        if(i==pd9)
        {for(int k=0;k<samples;k++)
        {
            {for(int j=0;j<46;j++)
            {freq9[(int)(telomother[k][j])]=freq9[(int)(telomother[k][j])]+1;}}
        }}
        
    }
        file<<endl;
        for(int pop=0;pop<(maxpd+1);pop++)
        {file<<population[pop]<<"\t";}
        
        
        
        
    }
    
    
    
    
    
    
    
    //results
    
    
   
    
    //histogram
    
    
    int hisgram=250;
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq1[j]=freq1[j]+freq1[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd1<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq1[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq2[j]=freq2[j]+freq2[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd2<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq2[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq3[j]=freq3[j]+freq3[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd3<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq3[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq4[j]=freq4[j]+freq4[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd4<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq4[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq5[j]=freq5[j]+freq5[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd5<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq5[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq6[j]=freq6[j]+freq6[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd6<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq6[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq7[j]=freq7[j]+freq7[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd7<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq7[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq8[j]=freq8[j]+freq8[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd8<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq8[j]<<"\t";
    }
    file<<endl;
    
    for(int j=1;j<lmax;j=j+hisgram)
    {
        for(int i=1;i<hisgram;i++)
        {freq9[j]=freq9[j]+freq9[j+i];}
    }
    
    file<<"pd"<<"\t"<<pd9<<endl;
    for(int j=1;j<16000;j=j+hisgram)
    {
        file<<freq9[j]<<"\t";
    }
    file<<endl;
    
    
    file.close();
    
    
    return 0;
}

