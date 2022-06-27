#include  <iostream>
#include  "stdio.h"
#include  <cmath>
#include  "getCommonSeq/ComSeq.h"
using namespace std;


void LoadM(const char *filename,sp_mat &V) {         
    vector<long long unsigned int> location_u;
    vector<long long unsigned int> location_m;
    vector<double> values;                    

    ifstream file(filename);                  
    int a, b;
    double c;                              
    while(file >> a >> b >> c) {                                   
        location_u.push_back(a);              
        location_m.push_back(b);              
        values.push_back(sqrt(c));                  
    }                                         

    umat lu(location_u);                      
    umat lm(location_m);                      
    umat location(join_rows(lm, lu).t());    
    sp_mat V1(location, vec(values)); 
    V = V1;                                        
}   

 






int main( int argc , char *argv[] )
{

   ComSeq test;
   char buffer[1024];

   FILE * ip =  fopen(argv[1],"r");
   FILE * op = fopen(argv[2],"w+");
   FILE * NameListIp=NULL;  
   FILE * configip=NULL;





   if( argc>2 )  
     NameListIp = fopen(argv[3],"r");


  

    for(int i =0;i<argc-1;i++)
    {
       if( strcmp(argv[i],"-config")==0)
       {
           configip = fopen(argv[i+1],"r");
       }
    }

    if( configip == NULL)
    {
        configip = fopen("config","r");
    }







   if((ip == NULL)||(op== NULL))
   {
      printf("usage: ./cluster MatrixFile  OuputFile  [ TableNameFile ]   [ -config  ConfigurationFile ]");
      printf("\n");
      return 0;
   }

    if( NameListIp != NULL)
    {
       test.ReadTableNameFile(NameListIp);
    }
   //test.matrixReadWithName(ip);
   test.LimitedRank = 1024;
   test.tolerance = 0.001;
   test.range = 0.12;
   test.rankSkip = true;
   test.Centering = false;
// test.matrixReadWithName(ip);
//   M_d.close();


   if(configip!= NULL)
   {


       while (!feof(configip))
       {
               fscanf(configip,"%[\b|\t]*",buffer);
               fscanf(configip,"%s",buffer);

               if(strcmp(buffer,"Tolerance")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.tolerance);
                   }
                   fscanf(ip,"%[\b|\t|\n]",buffer);
                   continue;
               }

            if(strcmp(buffer,"minRatio")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.minRatio);
                   }
                   fscanf(ip,"%[\b|\t|\n]",buffer);
                   continue;
               }




               if(strcmp(buffer,"ZeroTolerance")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.NoZeroRatio);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }

               if(strcmp(buffer,"VarianceThreshold")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.VarThreshold);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }






              if(strcmp(buffer,"CutOff")==0)
              {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.Log_Pvalue);
                      test.Log_Pvalue = log(test.Log_Pvalue);

                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }


                 if(strcmp(buffer,"MaxDepth")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%d",&test.MaxDepth);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }


                  if(strcmp(buffer,"MinSupport")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%d",&test.minsupport);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }

                    
                 if(strcmp(buffer,"CutThresholdRatio")==0)
               {
                   fscanf(configip,"%[\b|\t]*",buffer);
                   fscanf(configip,"%s",buffer);
                   test.CutThresHold = 0.0;

                   if(strcmp(buffer,"=")==0)
                   {
                      fscanf(configip,"%lf",&test.CutThresHold);
                   }
                   fscanf(configip,"%[\b|\t|\n]",buffer);
                   continue;
               }
       }

   }

   //test.matrixReadWithName(ip);
    sp_mat V;

   if(ip!=0)
   LoadM(argv[1],test.SampleS);


    test.SampleName.resize(test.SampleS.n_cols);

    if( NameListIp == NULL)
    {
         for(int j=0;j<test.SampleS.n_cols;j++)
        {
            test.SampleName[j]= std::to_string(j);
        }
    }
   
   if(test.SampleName.size()!=test.SampleS.n_cols)
   {
      printf("The length of table names is wrong.");
      return 0;
   }
   printf("Matrix have been read.\n");
   printf("n_rows: %d\n",test.SampleS.n_rows);
   printf("n_cols: %d\n",test.SampleS.n_cols);
   //test.filter();
   test.AttributeCluster();
   test.OutClusterWithNameCombact(op);
   return 0;

}
