#include<iostream.h>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
// Author Lalitha Viswanathan
// TCS (Tata Consultancy Services)
FILE* foutput=fopen("output.txt","w");
class DNA
{
	// Monomer: 4
	// Dimer: 4 x4 
	// etc.
	// hexa 4^(nmer size)
 int no_nuc[4],size,dimers[4][4],trimers[4][4][4],tetramers[4][4][4][4],pentamers[4][4][4][4][4],hexamers[4][4][4][4][4][4];
 int *nucleotide;
 int* DNA_seq;
 public:
	 // Constructor
	 DNA(char *filename);
	 // Destructor
	 ~DNA();
	 //sets counts for mon-hexa mers
	 void set_count();
	 void get_count(int) const;
	 //gets count for specific n-mer >6
	 void get_c(int);
	 int getsize() const{return size;}
	 int* getseq();
	 // Print DNA Seq
	 void DNA::printseq()
	 {
         for(int counter=0;counter<size;counter++)
	                 switch(DNA_seq[counter])
	                 {
                         case '0':
                             cout<<"A";
                             break;
                        case '1':
                             cout<<"U";
                            break;
                       case '2':
                            cout<<"G";
                            break;
                       case '3':
                           cout<<"C";
                           break;
               }
	 }
	 //copies sequence into a file
	 void DNA::copyseq(char *fname)
	 {
         FILE* fp;
         fp=fopen(fname,"w+");
         for(int counter=0;counter<size;counter++)
                 switch(DNA_seq[counter])
                {
                 case '0':
                      fputc('A',fp);
                      break;
                 case '1':
                      fputc('U',fp);
                       break;
                 case '2':
	              fputc('G',fp);
	              break;
                case '3':
                      fputc('C',fp);
                      break;
                 }
         fclose(fp);
	 }
	 //returns mono nucleotide count
 	int get_count(char nucleotide) const 
	{
	 if(nucleotide=='A')
		  return no_nuc[0];
	 else if(nucleotide=='T')
		  return no_nuc[1];
	 else if(nucleotide=='G')
		  return no_nuc[2];
	 else
		  return no_nuc[3];
 	}
	//given a di-hexamer string, returns count of that string in the DNA sequence
	int get_count(char* str) const
	{
		int index0,index1,index2,index3,index4,index5;
		if(strlen(str)>6)
			return -1;
		else
		{
			for(int counter=0;counter<strlen(str);counter++)
			{
					switch(str[counter])
					{
						case 'A':
						case 'a':
							if(counter==0)
								index0=0;
							else if(counter==1)
								index1=0;
							else if(counter==2)
								index2=0;
							else if(counter==3)
								index3=0;
							else if(counter==4)
								index4=0;
							else if(counter==5)
								index5=0;
							break;
						case 'T':
						case 't':
							if(counter==1)
								index0=1;
							else if(counter==1)
								index1=1;
							else if(counter==2)
								index2=1;
							else if(counter==3)
								index3=1;
							else if(counter==4)
								index4=1;
							else if(counter==5)
								index5=1;
							break;
						case 'G':
						case 'g':
							if(counter==2)
								index0=2;
							else if(counter==1)
								index1=2;
							else if(counter==2)
								index2=2;
							else if(counter==3)
								index3=2;
							else if(counter==4)
								index4=2;
							else if(counter==5)
								index5=2;
							break;
						case 'C':
						case 'c':
							if(counter==3)
								index0=3;
							else if(counter==1)
								index1=3;
							else if(counter==2)
								index2=3;
							else if(counter==3)
								index3=3;
							else if(counter==4)
								index4=3;
							else if(counter==5)
								index5=3;
							break;
					}
			}
			switch(strlen(str))
			{
				case 2:
					return dimers[index0][index1];
					break;
				case 3:
					return trimers[index0][index1][index2];
					break;
				case 4:
					return tetramers[index0][index1][index2][index3];
					break;
				case 5:
					return pentamers[index0][index1][index2][index3][index4];
					break;
				case 6:
					return hexamers[index0][index1][index2][index3][index4][index5];
					break;
			}
		}
	}
};

//function sets counts of mono,di,tri,tetra and hexamers 
void DNA::set_count()
{
       for(int i=0;i+5<size;i++)
       {
	 dimers[DNA_seq[i]][DNA_seq[i+1]]++;
         trimers[DNA_seq[i]][DNA_seq[i+1]][DNA_seq[i+2]]++;	 
         tetramers[DNA_seq[i]][DNA_seq[i+1]][DNA_seq[i+2]][DNA_seq[i+3]]++;	 
         pentamers[DNA_seq[i]][DNA_seq[i+1]][DNA_seq[i+2]][DNA_seq[i+3]][DNA_seq[i+4]]++;
         hexamers[DNA_seq[i]][DNA_seq[i+1]][DNA_seq[i+2]][DNA_seq[i+3]][DNA_seq[i+4]][DNA_seq[i+5]]++;	 
       }
return;       
}
//function sets count of n-mer where n>6
void DNA::get_c(int nmercount)
{
	nucleotide=(int *)malloc(sizeof(int)*(int)pow(4,nmercount));
			int num=0;
			cout<<"\nInside function size is "<<size<<"\n";
	for(int counter=0;counter+nmercount<size;counter++)
	{
			num=0;
		for(int innercounter=nmercount-1,k=counter;innercounter>=0;innercounter--,k++)
		{
			//num=0;
			num+=((int)pow(4,innercounter))*DNA_seq[k];
		}
	nucleotide[num]++;
	}
}
//fun for mapping integer to character
char map(int DNACounter)
{
	if(DNACounter==0)
		return 'A';
	else if(DNACounter==1)
		return 'T';
	else if(DNACounter==2)
		return 'G';
	else if(DNACounter==3)
		return 'C';
}
// Get NMer
char* get_nmer(int i,int desirednmerlength)
{
   char* a=(char *)malloc(sizeof(char)*desirednmerlength);
   for(int counterm=desirednmerlength-1,n=0;counterm>=0;counterm--,n++)
   {
		   a[n]=map(i/(int)pow(4,counterm));
		   i=i%(int)pow(4,counterm);
   }
  return a;
}
//returns count of n-mer
void DNA::get_count(int nmer) const
{
		
	if(nmer<=6)
	{
		switch(nmer)
		{
			case 1:
				for(int counter=0;counter<4;counter++)
						{
						cout<<map(counter)<<" : "<<no_nuc[counter]<<"\n";
						fprintf(foutput,"\n%c %d",map(counter),no_nuc[counter]);					
						}
				break;
			case 2:
				for(int counter=0;counter<4;counter++)
				for(int innercounter=0;innercounter<4;innercounter++)
				{
						cout<<map(counter)<<map(innercounter)<<" : "<<dimers[counter][innercounter]<<"\n";
					fprintf(foutput,"\n%c%c %d",map(counter),map(innercounter),dimers[counter][innercounter]);					
				}
				break;
			case 3:
				for(int counter=0;counter<4;counter++)
				for(int innercounter=0;innercounter<4;innercounter++)
				for(int kcounter=0;kcounter<4;kcounter++)
				{
						cout<<map(counter)<<map(innercounter)<<map(kcounter)<<" : "<<trimers[counter][innercounter][kcounter]<<"\n";
						fprintf(foutput,"\n%c%c%c %d",map(counter),map(innercounter),map(kcounter),trimers[counter][innercounter][kcounter]);					
				}
				break;
			case 4:
				for(int i=0;i<4;i++)
				for(int j=0;j<4;j++)
				for(int k=0;k<4;k++)
				for(int l=0;l<4;l++)
				{
						cout<<map(i)<<map(j)<<map(k)<<map(l)<<" : "<<tetramers[i][j][k][l]<<"\n";
						fprintf(foutput,"\n%c%c%c%c %d",map(i),map(j),map(k),map(l),tetramers[i][j][k][l]);						  }
				break;
			case 5:
				for(int counteri=0;counteri<4;counteri++)
				for(int counterj=0;counterj<4;counterj++)
				for(int counterk=0;counterk<4;counterk++)
				for(int counterl=0;counterl<4;counterl++)
				for(int counterm=0;counterm<4;counterm++)
				{
					cout<<map(counteri)<<map(counterj)<<map(counterk)<<map(counterl)<<map(counterm)<<" : "<<pentamers[counteri][counterj][counterk][counterl][counterm]<<"\n";
					fprintf(foutput,"\n%c%c%c%c%c %d",map(counteri),map(counterj),map(counterk),map(counterl),map(counterm),pentamers[counteri][counterj][counterk][counterl][counterm]);						  
				}
				break;
			case 6:
				for(int counteri=0;counteri<4;counteri++)
				for(int counterj=0;counterj<4;counterj++)
				for(int counterk=0;counterk<4;counterk++)
				for(int counterl=0;counterl<4;counterl++)
				for(int counterm=0;counterm<4;counterm++)
				for(int countern=0;countern<4;countern++)
				{
					cout<<map(counteri)<<map(counterj)<<map(counterk)<<map(counterl)<<map(counterm)<<map(countern)<<" : "<<hexamers[counteri][counterj][counterk][counterl][counterm][countern]<<"\n";
					fprintf(foutput,"\n%c%c%c%c%c%c %d",map(counteri),map(counterj),map(counterk),map(counterl),map(counterm),map(countern),hexamers[counteri][counterj][counterk][counterl][counterm][countern]);						  }
				break;
		}
			
	}
	else
	{
		for(int counteri=0;counteri<pow(4,nmer);counteri++)
		{
				char* z=get_nmer(counteri,nmer);
				cout<<z <<" : "<<nucleotide[counteri]<<"\n";
				fprintf(foutput,"\n%s %d",z,nucleotide[counteri]);
		}
	}
}
//Protein class
class protein
{
 int protein_count[22];	
 char* protein_seq[6];
 int size[6];
 public:
  protein(DNA);
  void get_count()const;
  void protein::printseq(int );
  void protein::copyseq(char*);
  int get_count(char str)const;
};

// Constructor
 protein::protein(DNA DNAclass)
 {
	 int counter,innercounter,*sequence;
	 sequence=DNAclass.getseq();
	 int tempsize=DNAclass.getsize(); 
	 for(innercounter=0;innercounter<3;innercounter++)
	 {
	  protein_seq[innercounter]=(char*)malloc(sizeof(char)*(tempsize/3));
	  size[innercounter]=0;
	  for(counter=innercounter;(counter+3)<tempsize;counter+=3)
  	  {
		  int number=(sequence[counter]*100)+(sequence[counter+1]*10)+sequence[counter+2];
		  switch(number)
		  {
			  case 111:
			  case 113:
				  protein_seq[innercounter][size[innercounter]++]='F';
				  protein_count[0]++;
				  break;
			  case 110:
			  case 112:
			  case 311:
			  case 313:
			  case 310:
			  case 312:
			          protein_seq[innercounter][size[innercounter]++]='L';
				  protein_count[1]++;
				  break;
			  case 131:
			  case 133:
			  case 130:
			  case 132:
			  case 021:		  
			  case 023:		  
			          protein_seq[innercounter][size[innercounter]++]='S';
				  protein_count[2]++;
				  break;
			  case 101:
			  case 103:
				  protein_seq[innercounter][size[innercounter]++]='Y';
				  protein_count[3]++;
				  break;
			  case 100:
			  case 102:
			  case 120:
			          protein_seq[innercounter][size[innercounter]++]='X';
				  protein_count[4]++;
			          break;
			  case 121:
			  case 123:
			          protein_seq[innercounter][size[innercounter]++]='C';
				  protein_count[5]++;
			          break;
			  case 122:
			           protein_seq[innercounter][size[innercounter]++]='W';
				   protein_count[6]++;
			           break;
			  case 331:
			  case 333:
			  case 330:
			  case 332:
				   protein_seq[innercounter][size[innercounter]++]='P';
				   protein_count[7]++;
				   break;
			  case 301:
			  case 303:
				   protein_seq[innercounter][size[innercounter]++]='H';
				   protein_count[8]++;
				   break;
			  case 300:
			  case 302:
				   protein_seq[innercounter][size[innercounter]++]='Q';
				   protein_count[9]++;
				   break;
			  case 321:
			  case 323:
			  case 320:
			  case 322:
			  case 20:		  
			  case 22:
				   protein_seq[innercounter][size[innercounter]++]='R';
				   protein_count[10]++;
				   break;
		          		   
			  case 11:		  
			  case 13:		  
			  case 10:
		 		   protein_seq[innercounter][size[innercounter]++]='I';
				   protein_count[11]++;
				   break;		   
			  case 12:
		 		   protein_seq[innercounter][size[innercounter]++]='M';
				   protein_count[12]++;
				   break;		   
			  case 31:		  
			  case 33:		  
			  case 30:		  
			  case 32:
		 		  protein_seq[innercounter][size[innercounter]++]='T';
				   protein_count[13]++;
       				  break;				  
			  case 1:		  
			  case 3:
				   protein_seq[innercounter][size[innercounter]++]='N';
				   protein_count[14]++;
		 		   break;		   
			  case 0:		  
			  case 2:
		 		   protein_seq[innercounter][size[innercounter]++]='K';
				   protein_count[15]++;
				   break;		   
		          case 211:		   
		          case 213:		   
		          case 210:		   
		          case 212:
				   protein_seq[innercounter][size[innercounter]++]='V';
				   protein_count[16]++;
				   break;		   
		          case 231:		   
		          case 233:		   
		          case 230:		   
		          case 232:
				   protein_seq[innercounter][size[innercounter]++]='A';
				   protein_count[17]++;
				   break;		   
		          case 201:		   
		          case 203:
				   protein_seq[innercounter][size[innercounter]++]='D';
				   protein_count[18]++;
				   break;		   
		          case 200:		   
		          case 202:
				   protein_seq[innercounter][size[innercounter]++]='E';
				   protein_count[19]++;
				   break;			   
		          case 221:		   
		          case 223:		   
		          case 220:		   
		          case 222:
				   protein_seq[innercounter][size[innercounter]++]='G';
				   protein_count[20]++;
				   break;			   
		  }
	  
	  }
	 cout<<"After generating protein sequence for the forward strand in 3 reading frames:"<<innercounter<<endl;
	 }
	 for(innercounter=3;innercounter<6;innercounter++)
	 {
	  protein_seq[innercounter]=(char*)malloc(sizeof(char)*(DNAclass.getsize()/3));
	  size[innercounter]=0;
	  for(int i=DNAclass.getsize()-innercounter+2;i-3>0;i-=3)
  	  {
		  int number=(sequence[i]*100)+(sequence[i+1]*10)+sequence[i+2];
		  switch(number)
		  {
			  case 111:
			  case 311:
				  protein_seq[innercounter][size[innercounter]++]='F';
				   protein_count[0]++;
				  break;
			  case 011:
			  case 211:
			  case 113:
			  case 313:
			  case 013:
			  case 213:
			          protein_seq[innercounter][size[innercounter]++]='L';
				   protein_count[1]++;
				  break;
			  case 131:
			  case 331:
			  case 031:
			  case 231:
			  case 120:		  
			  case 320:		  
			          protein_seq[innercounter][size[innercounter]++]='S';
				   protein_count[2]++;
				  break;
			  case 101:
			  case 301:
				  protein_seq[innercounter][size[innercounter]++]='Y';
				   protein_count[3]++;
				  break;
			  case 001:
			  case 201:
			  case 021:
			          protein_seq[innercounter][size[innercounter]++]='X';
				   protein_count[4]++;
			          break;
			  case 121:
			  case 321:
			          protein_seq[innercounter][size[innercounter]++]='C';
				   protein_count[5]++;
			          break;
			  case 221:
			           protein_seq[innercounter][size[innercounter]++]='W';
				   protein_count[6]++;
			           break;
			  case 133:
			  case 333:
			  case 033:
			  case 233:
				   protein_seq[innercounter][size[innercounter]++]='P';
				   protein_count[7]++;
				   break;
			  case 103:
			  case 303:
				   protein_seq[innercounter][size[innercounter]++]='H';
				   protein_count[8]++;
				   break;
			  case 003:
			  case 203:
				   protein_seq[innercounter][size[innercounter]++]='Q';
				   protein_count[9]++;
				   break;
			  case 123:
			  case 323:
			  case 023:
			  case 223:
			  case 020:		  
			  case 220:
				   protein_seq[innercounter][size[innercounter]++]='R';
				   protein_count[10]++;
				   break;
		          		   
			  case 110:		  
			  case 310:		  
			  case 010:
		 		   protein_seq[innercounter][size[innercounter]++]='I';
				   protein_count[11]++;
				   break;		   
			  case 210:
		 		   protein_seq[innercounter][size[innercounter]++]='M';
				   protein_count[12]++;
				   break;		   
			  case 130:		  
			  case 330:		  
			  case 030:		  
			  case 230:
		 		  protein_seq[innercounter][size[innercounter]++]='T';
				   protein_count[13]++;
       				  break;				  
			  case 100:		  
			  case 300:
				   protein_seq[innercounter][size[innercounter]++]='N';
				   protein_count[14]++;
		 		   break;		   
			  case 000:		  
			  case 200:
		 		   protein_seq[innercounter][size[innercounter]++]='K';
				   protein_count[15]++;
				   break;		   
		          case 112:		   
		          case 312:		   
		          case 012:		   
		          case 212:
				   protein_seq[innercounter][size[innercounter]++]='V';
				   protein_count[16]++;
				   break;		   
		          case 132:		   
		          case 332:		   
		          case 032:		   
		          case 232:
				   protein_seq[innercounter][size[innercounter]++]='A';
				   protein_count[17]++;
				   break;		   
		          case 102:		   
		          case 302:
				   protein_seq[innercounter][size[innercounter]++]='D';
				   protein_count[18]++;
				   break;		   
		          case 002:		   
		          case 202:
				   protein_seq[innercounter][size[innercounter]++]='E';
				   protein_count[19]++;
				   break;			   
		          case 122:		   
		          case 322:		   
		          case 022:		   
		          case 222:
				   protein_seq[innercounter][size[innercounter]++]='G';
				   protein_count[20]++;
				   break;			   
		  }
	  
	  }
	 cout<<"After generating protein sequence for the reverse strand in 3 reading frames:"<<innercounter<<endl;
	 }
 }
//returns count of all amino acids
void protein::get_count() const
{
	for(int counter=0;counter<20;counter++)
	{
		switch(counter)
		{
			case 0:
				cout<<"F :"<<protein_count[0]<<endl;
				fprintf(foutput,"\nF : %d",protein_count[0]);
				break;
			case 1:
				cout<<"L :"<<protein_count[1]<<endl;
				fprintf(foutput,"\nL : %d",protein_count[1]);
				break;
			case 2:
				cout<<"S :"<<protein_count[2]<<endl;
				fprintf(foutput,"\nS : %d",protein_count[2]);
				break;
			case 3:
				cout<<"Y :"<<protein_count[3]<<endl;
				fprintf(foutput,"\nY : %d",protein_count[3]);
				break;
			case 4:
				cout<<"X :"<<protein_count[4]<<endl;
				fprintf(foutput,"\nX : %d",protein_count[4]);
				break;
			case 5:
				cout<<"C :"<<protein_count[5]<<endl;
				fprintf(foutput,"\nC : %d",protein_count[5]);
				break;
			case 6:
				cout<<"W :"<<protein_count[6]<<endl;
				fprintf(foutput,"\nW : %d",protein_count[6]);
				break;
			case 7:
				cout<<"P :"<<protein_count[7]<<endl;
				fprintf(foutput,"\nP : %d",protein_count[7]);
				break;
			case 8:
				cout<<"H :"<<protein_count[8]<<endl;
				fprintf(foutput,"\nH : %d",protein_count[8]);
				break;
			case 9:
				cout<<"Q :"<<protein_count[9]<<endl;
				fprintf(foutput,"\nQ : %d",protein_count[9]);
				break;
			case 10:
				cout<<"R :"<<protein_count[10]<<endl;
				fprintf(foutput,"\nR : %d",protein_count[10]);
				break;
			case 11:
				cout<<"I :"<<protein_count[11]<<endl;
				fprintf(foutput,"\nI : %d",protein_count[11]);
				break;
			case 12:
				cout<<"M :"<<protein_count[12]<<endl;
				fprintf(foutput,"\nM : %d",protein_count[12]);
				break;
			case 13:
				cout<<"T :"<<protein_count[13]<<endl;
				fprintf(foutput,"\nT : %d",protein_count[13]);
				break;
			case 14:
				cout<<"N :"<<protein_count[14]<<endl;
				fprintf(foutput,"\nN : %d",protein_count[14]);
				break;
			case 15:
				cout<<"K :"<<protein_count[15]<<endl;
				fprintf(foutput,"\nK : %d",protein_count[15]);
				break;
			case 16:
				cout<<"V :"<<protein_count[16]<<endl;
				fprintf(foutput,"\nV : %d",protein_count[16]);
				break;
			case 17:
				cout<<"A :"<<protein_count[17]<<endl;
				fprintf(foutput,"\nA : %d",protein_count[17]);
				break;
			case 18:
				cout<<"D :"<<protein_count[18]<<endl;
				fprintf(foutput,"\nD : %d",protein_count[18]);
				break;
			case 19:
				cout<<"E :"<<protein_count[19]<<endl;
				fprintf(foutput,"\nE : %d",protein_count[19]);
				break;
			case 20:
				cout<<"G :"<<protein_count[20]<<endl;
				fprintf(foutput,"\nG : %d",protein_count[20]);
				break;
		}
	}
}
//returns count of particular amino acid
int protein::get_count(char AA)const
{
	switch(AA)
	{
		case 'F':
			fprintf(foutput,"\nF : %d",protein_count[0]);
			return protein_count[0];
			break;
		case 'L':
			fprintf(foutput,"\nL : %d",protein_count[1]);
			return protein_count[1];
			break;
		case 'S':
			fprintf(foutput,"\nS : %d",protein_count[2]);
			return protein_count[2];
			break;
		case 'Y':
			fprintf(foutput,"\nY : %d",protein_count[3]);
			return protein_count[3];
			break;
		case 'X':
			fprintf(foutput,"\nX : %d",protein_count[4]);
			return protein_count[4];
			break;
		case 'C':
			fprintf(foutput,"\nC : %d",protein_count[5]);
			return protein_count[5];
			break;
		case 'W':
			fprintf(foutput,"\nW : %d",protein_count[6]);
			return protein_count[6];
			break;
		case 'P':
			fprintf(foutput,"\nP : %d",protein_count[7]);
			return protein_count[7];
			break;
		case 'H':
			fprintf(foutput,"\nH : %d",protein_count[8]);
			return protein_count[8];
			break;
		case 'Q':
			fprintf(foutput,"\nQ : %d",protein_count[9]);
			return protein_count[9];
			break;
		case 'R':
			fprintf(foutput,"\nR : %d",protein_count[10]);
			return protein_count[10];
			break;
		case 'I':
			fprintf(foutput,"\nI : %d",protein_count[11]);
			return protein_count[11];
			break;
		case 'M':
			fprintf(foutput,"\nM : %d",protein_count[12]);
			return protein_count[12];
			break;
		case 'T':
			fprintf(foutput,"\nT : %d",protein_count[13]);
			return protein_count[13];
			break;
		case 'N':
			fprintf(foutput,"\nN : %d",protein_count[14]);
			return protein_count[14];
			break;
		case 'V':
			fprintf(foutput,"\nV : %d",protein_count[15]);
			return protein_count[15];
			break;
		case 'A':
			fprintf(foutput,"\nA : %d",protein_count[16]);
			return protein_count[16];
			break;
		case 'D':
			fprintf(foutput,"\nD : %d",protein_count[17]);
			return protein_count[17];
			break;
		case 'E':
			fprintf(foutput,"\nE : %d",protein_count[18]);
			return protein_count[18];
			break;
		case 'G':
			fprintf(foutput,"\nG : %d",protein_count[19]);
			return protein_count[19];
			break;
	}
}
// DNA Class Destructor
 DNA::~DNA()
 {
 }
//prints protein sequence
void protein::printseq(int lengthofproteinseq)
{
	        for(int i=0;i<size[lengthofproteinseq];i++)
	                cout<<protein_seq[i];
}
//copies protein sequence into a file
void protein::copyseq(char *fname)
{
	        FILE* fp;
	        fp=fopen(fname,"w+");
	        for(int counterj=0;counterj<6;counterj++)
		        for(int counteri=0;counteri<size[counterj];counteri++)
		                fputc(protein_seq[counterj][counteri],fp);
		        fclose(fp);
}
//returns DNA sequence
int* DNA::getseq()
{
	return DNA_seq;
}
//DNA constructor
 DNA::DNA(char* filename)
 {
 	  FILE* fp;
	  char str[500],ch;
	  int temp=0;
	 if( !(fp=fopen(filename ,"r")))
	 {
			 printf("Filename %s does not exist\n",filename);
			 exit(1);
	 }
	  fseek(fp,0L,SEEK_END);
	  size=ftell(fp);
	  DNA_seq=(int *)(malloc(sizeof(int)*size));
	  fseek(fp,0L,SEEK_SET);
	  fgets(str,500,fp);
	 while(!feof(fp))
	  {	
	     ch=fgetc(fp);
	     if((ch=='A')||(ch=='a'))
	     	    {
		    	no_nuc[0]++;
   	                DNA_seq[temp++]=0;
		    }
	    else if((ch=='T')||(ch=='t'))
	    	    {
		    	no_nuc[1]++;
	            	DNA_seq[temp++]=1;
		    }
	    else if((ch=='G')||(ch=='g'))
	    	    {
		    	no_nuc[2]++;
	            	DNA_seq[temp++]=2;
		    }
	    else if((ch=='C')||(ch=='c'))
	    	    {
			no_nuc[3]++;
	                DNA_seq[temp++]=3;
		    }
  	  }
          fclose(fp);
          cout<<"Size of the genome ="<<size<<"bp\n";
 }

 // Main program 
main(int argc, char* argv[])
{
       	if(argc<2)
	 {
	  printf("enter the input in the form <filename1>....<filename n>\n"); 
	  exit(1);
 	 }
		//finds composition for n number of files
	  for(int y=1;y<argc;y++)
	  {
	  int count,choice=1,beg_range,end_range;
      DNA DNAsequence(argv[y]);
	  //sets counts for mono-hexa mers for DNA sequence
	  DNAsequence.set_count();
	  cout<<"\nEnter the n-mer count u wish to observe: ";
	 		  cin>> count;
			  if(count>6)
			  {
			  DNAsequence.get_c(count);
			  DNAsequence.get_count(count);
			  }
			  else
					  DNAsequence.get_count(count);
          protein PROTsequence(DNAsequence);
	  cout<<"\n Press 1 to view count of all amino acids in the sequence or press 2 to view count of particular amino acid";
	  int choice2;
	  cin>>choice2;
	  if((choice2<1) ||(choice2>2))
		  choice2=1;
	  switch(choice2)
	  {
		  case 1:
			  PROTsequence.get_count();
			  break;
		  case 2:
			  char aa;
			  cin>>aa;
			  PROTsequence.get_count(aa);
			  break;
	  }
	  }
fclose(foutput);
}
