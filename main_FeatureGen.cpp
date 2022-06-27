#include  <iostream>
#include  "stdio.h"
#include  <cmath>
#include  "getCommonSeq/ComSeq.h"
using namespace std;



int main( int argc , char *argv[] )
{
          
     ClusterToFeature Test;



     FILE * FileName=NULL; 


     for(int i =0;i<argc-1;i++)
     {
       if( strcmp(argv[i],"-N")==0)
       {
           FileName = fopen(argv[i+1],"r");
       }
     }

          
     if(argv[1]==NULL)
     {  
          printf("FeatureGen ToFeature\n");
          printf("FeatureGen Decipher\n");

          return 0;
      }

     if( strcmp(argv[1],"ToFeature")==0)
     {


        FILE * ip =  fopen(argv[2],"r");
        FILE * op =  fopen(argv[3],"w+");

        if(ip==NULL||op==NULL||argc<5)
        {
            printf("FeatureGen ToFeature   InputFile OutputFIle AverageCorrelation AverageComponentRatio\n");
            return 0;
        }
        Test.tolerance = 0.01;
        Test.CorThreshold = atof(argv[4]);
        Test.CompRThreshold = atof(argv[5]);
        Test.ClusterRead(ip);
        Test.OutPutFeature(op);
     }

     if( strcmp(argv[1],"Decipher")==0)
     {


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op =  fopen(argv[4],"w+");

         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen Decipher FeatureFile   InputFile -N TableNameFile OutputFile\n");
            return 0;
         }

        Test.ReadFeature(ip);
        fclose(ip1);
        Test.SampleReadS(argv[3]);


        if(FileName==NULL)
        {
             Test.SampleName.resize(Test.SampleS.n_cols);

            for(int j=0;j<Test.SampleS.n_cols;j++)
            {
                Test.SampleName[j]= std::to_string(j);
            }
        }else
         {
                Test.ReadSNameFile(FileName);
         }

         if(Test.SampleName.size()!=Test.SampleS.n_cols)
         {
             printf("The length of sample name list is wrong!\n");
             return(0);
         }

        Test.SampleDecipherS(op);
     }



 }
