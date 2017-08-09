#include<stdlib.h>
#include<stdio.h>
//#include<conio.h>
#include<time.h>
#include<cstring>
#include<math.h>
#include <iostream>

#define pi 3.14159
using namespace std;

/*system parameters*/
const int S=5;    //station nubmber
const int B=60;    //buses number
int C=80;    //capacity of each bus
int C0=40;  //seat capacity
double tbar=3;	//the mean time of travel between stations
time_t  start;
time_t  sample;

/*Energy parameters*/

//fuel: is modelede like k*x, in which k is somehow a stiffness and x in the excess velocity realtive to mean volicity
double x=1;	//difference between velocity and mean velocity
double K=50;	//stiffness
//comfort: is modeled as a two value parameter. 1 if a person has a seat and 0 if he/she hasn't. 
//The total comfort (which is stored in "comf" is average of comfort for all passengers
double comf=0;  //comfort variable


/*specific time parameters*/
double tin=0.05;    //time needed for each person to get on
double tout=0.05;  //time needed for each person to get off
double tb=1;    //minimum Headway time
double ts=1;    //minimum time from departure (from a station) to arriving (to another station)
double tline=5; //lineman dt (the time interval by which the lineman enters a bus into the line)
double L0=10;    //landa scale
double T=4.5*B*ts;       //work time:landa will be 0 after T
double kTmax=1000;   //maximum Temperature
double kT;	//Temperature
int imax=100;  //counter max for evaluating r

/*scaler variables*/
double r;       //acceptance ratio
double R;       //cumulative r
double Ei;      //Initial Energy
double Ef;      //Final Energy
double E;       //Energy
double dE;      //delta Energy
double Emin=100000000;	//It is a dummy variable. we store the minimum Energy throughout our survey in phase space in this variable.
double q=.5;	//This guy is the one who changes t elements.

/*variables*/
double zeroarray[B+2][S+2]={0};
double Nin[B+2][S+2]={0} ;  //Number of passengers in bus
double n[B+2][S+2]={0} ;  //number of passengers in station: in the moment of bus's arrival
double nn[B+2][S+2]={0} ; //number of passengers in station: after leaving the bus
double t[B+2][S+2]={0};   //schedule
double tt[B+2][S+2]={0} ; //departure schedule
double En[B+2][S+2]={0} ;  //Energy matrix

/*copies*/
double Nini[B+2][S+2];   //Number of passengers in bus
double ni[B+2][S+2];   //number of passengers in station: in the moment of bus's arrival
double nni[B+2][S+2];  //number of passengers in station: after leaving the bus
double ti[B+2][S+2];   //schedule
double tti[B+2][S+2];  //departure schedule
double Eni[B+2][S+2];     //energy matrix
double landai[B+2][S+2];	//landa matrix

/*inputs*/
double landa[B+2][S+2];   //rate of entering the station
double beta[S+2];    //probability of getting off of the bus

/*spare variables*/
double cond;	//condition: This evaluates 
double delta=T/5; //variation in time (as a guess it is set to be sth)



/*LANDA function*/
double Landa (int station,int bus, double t){
    if (t<T){
        double tstar=(T-10*t)/T;
        landa[bus][station]=L0*(sin(pi/(S)*station))*(sin(pi/(T)*t));//(pow(3*t/T,2))*exp(-pow(tstar,2));
        }
    else{
        landa[bus][station]=0;
    }
    //cout<<"b0,s0,landa "<<bus<<"\t"<<station<<"\t"<<(landa[bus][station])<<endl;
}

/*BETA function*/
double Beta(){
    for (int s=1;s<S+1;s++){
        beta[s]=sin(pi/2*double(s)/S);
        //cout<<beta[s]<<endl;
    }
}


/*TDET function*/
int Tdet(int k){
    int Shart=1;	//This variable checks if a bus with number "b-k" is out of line or not.

    for (int b=1; b<B+1;b++){
        for (int s=1; s<S+1;s++){

            //setting the first element the beginning of the time
            if (s==1 ){
                t[b][s]=t[b-1][s]+tline;
                if (b>k){
                    Shart=t[b][s]>tt[b-k][S];
                    while (Shart==0){
                        t[b][s]+= (q)/4;      //it just increases t
                        Shart=t[b][s]>tt[b-k][S];
                        }
                    }
                }
            else{
                t[b][s]= max(tt[b-1][s]+10*tb,tt[b][s-1]+10*ts);    //This is just an arbitrary t[b][s] which also does satisfy constraints
                //cout<<t[b][s]<<endl;
                }
            

            Landa(s,b,t[b][s]);     //computes landa of s and b
            //cout<<s<<"\t"<<b<<endl;
            n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
            cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
            if (cond>=n[b][s]){
                    nn[b][s]=0;
                    Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                    }
            else{
                nn[b][s]=n[b][s]-cond;
                Nin[b][s+1]=C;
                }

            tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
            //cout<<t[b][s]<<"\t"<<tt[b][s]<<"\t"<<tin<<"\t"<<tout<<endl;
            //cout<<Nin[b][s]<<"\t"<<beta[s]<<"\t"<<endl;
            }
        }
    
}

/*ENERGY function*/
int Energy(int b0,int s0,int k,double tbar){
    //cout<<"In ENERGY b0 is "<<b0<<" and s0 is "<<s0<<endl;
	
	
	//evaluating Energy for bus b0
    for (int s=s0;s<S+1;s++){
        
		//evaluating comfort
        if (Nin[b0][s]<=C0){
            comf=1;
            }
        else{
            comf=double(C0)/Nin[b0][s];
            }
            
		//evaluating velocity excess
		x=t[b0][s]-t[b0][s-1]-tbar;
		if (x>=0){
			x=0;	//no consumption/erosion when speed is not high.
		}
		
        if (s!=1){
        	
            En[b0][s]= En[b0][s-1] + nn[b0-1][s]*(t[b0][s]-t[b0-1][s]) + landa[b0][s]/2*(t[b0][s]-t[b0-1][s])*(t[b0][s]-t[b0-1][s]) + Nin[b0][s]*(t[b0][s]-t[b0][s-1])/comf*(1-K*x);
            }
        else{
            En[b0][s]= En[b0-1][S] + nn[b0-1][s]*(t[b0][s]-t[b0-1][s]) + landa[b0][s]/2*(t[b0][s]-t[b0-1][s])*(t[b0][s]-t[b0-1][s]) + Nin[b0][s]*(t[b0][s]-t[b0][s-1])/comf*(1-K*x);
            }
        //cout<<"in ENERGY k is"<<k<<" E of "<<b0<<","<<s<<"\t"<<En[b0][s]<<endl;
        }

    //evaluating Energy for buses with number greater than b0
    for (int b=b0+1;b<B+1;b++){
        for (int s=1;s<S+1;s++){
            //evaluating comfort
            if (Nin[b][s]<=C0){
                comf=1;
                }
            else{
                comf=double(C0)/Nin[b][s];
                }
			
			//evaluating velocity excess
			x=t[b0][s]-t[b0][s-1]-tbar;
			if (x>=0){
				x=0;	//no consumption/erosion when speed is not high.
				}
		
            if (s!=1){
                En[b][s]= En[b][s-1] + nn[b-1][s]*(t[b][s]-t[b-1][s]) + landa[b][s]/2*(t[b][s]-t[b-1][s])*(t[b][s]-t[b-1][s]) + Nin[b][s]*(t[b][s]-t[b][s-1])/comf*(1-K*x);
                }
            else{
                En[b][s]= En[b-1][S] + nn[b-1][s]*(t[b][s]-t[b-1][s]) + landa[b][s]/2*(t[b][s]-t[b-1][s])*(t[b][s]-t[b-1][s]) + Nin[b][s]*(t[b][s]-t[b][s-1])/comf*(1-K*x);
                }
            //cout<<"in ENERGY k is"<<k<<" E of "<<b<<","<<s<<"\t"<<En[b][s]<<endl;
            }
		}
    

}


// FUNC function
int Func(int b0,int s0,int k){
    int Shart=1;
	
	
	//cout<<"in FUNC no change has been made. b and s are\t"<<b0<<"\t"<<s0<<"\ttt and t are "<<tt[b0][s0]<<"\t"<<t[b0][s0]<<endl;
	//Variation in t[b0][s0]:...
    t[b0][s0]+= 2*delta*(((double)(rand()/(double)RAND_MAX))-0.5);
	
    //checking SHART= just k trains are in the line
    if (b0>k && s0==1){
        Shart= t[b0][1]>tt[b0-k][S];
        }
        
    //check whether the constraints are satisfied or not
    while (t[b0][s0]-tt[b0-1][s0]<tb || t[b0][s0]-tt[b0][s0-1]<ts || Shart==0 ){
        t[b0][s0]+= q;      //it just increases t (no decreasing as it is not large enough!) 
		
        //checking SHART= just k trains are in the line
        if (b0>k && s0==1){
            Shart= (t[b0][s0]>tt[b0-k][S]); 
            }
        }
    
	
    Landa(s0,b0,t[b0][s0]);     //computes landa of s and b
    n[b0][s0]=nn[b0-1][s0]+landa[b0][s0]*(t[b0][s0]-tt[b0-1][s0]);   //people on station when b arrives s
    cond=C-Nin[b0][s0]*(1-beta[s0]);   //capacity of the bus
    if (cond>=n[b0][s0]){
        nn[b0][s0]=0;
        Nin[b0][s0+1]=Nin[b0][s0]*(1-beta[s0])+n[b0][s0];
        }
    else{
        nn[b0][s0]=n[b0][s0]-cond;
        Nin[b0][s0+1]=C;
        }
    tt[b0][s0]=t[b0][s0]+tout*beta[s0]*Nin[b0][s0]+tin*(Nin[b0][s0+1]-(1-beta[s0])*Nin[b0][s0]);  //time of departure
	//cout<<"in FUNC up to here, just the variation has been made. b and s are\t"<<b0<<"\t"<<s0<<"\ttt and t are "<<tt[b0][s0]<<"\t"<<t[b0][s0]<<endl;
	//cout<<"in FUNC b and s and k "<<b0<<" "<<s0<<" "<<k<<"\tt[b][s]="<<t[b0][s0]<<"\ttt[b-k][S]="<<tt[b0-k][S]<<endl;
	
	
	
    //Variation in the only block that is influenced by t[b0][s0]
    if (b0+k>B){
		
        for (int b=b0;b<B+1;b++){
            if (b==b0){
                for (int s=s0+1; s<S+1;s++){
					
                    t[b][s]=max(tt[b-1][s]+tb,tt[b][s-1]+ts);
					//printf("line=238 b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
                    //checking SHART= just k trains are in the line
                    if (b>k && s==1){
                        Shart= t[b][s]>tt[b-k][S];
                        
                        }
                    //check whether the constraints are satisfied of not
                    while (Shart==0){
                        t[b][s]+= q;     //it just increases t

                        //checking SHART= just k trains are in the line
                        if (b>k && s==1){
                            Shart= t[b][s]>tt[b-k][S];
                            }
                        }

                    Landa(s,b,t[b][s]);     //computes landa of s and b
                    n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
                    cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
                    if (cond>=n[b][s]){
                        nn[b][s]=0;
                        Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                        }
                    else{
                        nn[b][s]=n[b][s]-cond;
                        Nin[b][s+1]=C;
                        }
                    tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
                    }
                }
            else{
                for (int s=s0; s<S+1;s++){

                    t[b][s]=max(tt[b-1][s]+tb,tt[b][s-1]+ts);
					//printf("line=272 b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
                	
                    //checking SHART= just k trains are in the line
                    if (b>k && s==1){
                        Shart= t[b][s]>tt[b-k][S];
                        }
                    //check whether the constraints are satisfied of not
                    while (Shart==0 ){
                        t[b][s]+= q;      //it just increases t 
                        
                        //checking SHART= just k trains are in the line
                        if (b>k && s==1){
                            Shart= t[b][s]>tt[b-k][S];
                            }
                        }
                        

                    Landa(s,b,t[b][s]);     //computes landa of s and b
                    n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
                    cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
                    if (cond>=n[b][s]){
                        nn[b][s]=0;
                        Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                        }
                    else{
                        nn[b][s]=n[b][s]-cond;
                        Nin[b][s+1]=C;
                        }
                    tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
                    }
            	}
        	}
    	}

    //Variation in the two block  that are influenced by t[b0][s0]: (b0+k<B)
    else{
		
        //The Square Block
        for (int b=b0; b<b0+k;b++){
            if (b==b0){
                for (int s=s0+1; s<S+1;s++){

                    t[b][s]=max(tt[b-1][s]+tb,tt[b][s-1]+ts);
					//printf("line=315 b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
                    //checking SHART= just k trains are in the line
                    if (b>k && s==1){
                        Shart= t[b][s]>tt[b-k][S];
                        }
                    //check whether the constraints are satisfied of not
                    while (Shart==0 ){
                        t[b][s]+= q;//*delta/2;      //it just increases t (no decreasing as it is not large enough!) but the interval is same as the positive interval of randomizing t

                        //checking SHART= just k trains are in the line
                        if (b>k && s==1){
                            Shart= t[b][s]>tt[b-k][S];
                            }
                        }

                    Landa(s,b,t[b][s]);     //computes landa of s and b
                    n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
                    cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
                    if (cond>=n[b][s]){
                        nn[b][s]=0;
                        Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                        }
                    else{
                        nn[b][s]=n[b][s]-cond;
                        Nin[b][s+1]=C;
                        }
                    tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
                    }  
                }
            else{
            	//cout<<k<<endl;
                for (int s=s0; s<S+1;s++){

                    t[b][s]=max(tt[b-1][s]+tb,tt[b][s-1]+ts);
					//printf("b= %d and s=%d and k=%d and t[b][s]=%f and tt[b][s]=%f and tt[b-k][S]=%f\n",b,s,k,t[b][s],tt[b][s],tt[b-k][S]);
                    //printf("line=350 b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
					
					//checking SHART= just k trains are in the line
                    if (b>k && s==1){
                        Shart= t[b][s]>tt[b-k][S];
                        }
                    //check whether the constraints are satisfied of not
                    while (Shart==0 ){
                        t[b][s]+= q;//*delta/2;      //it just increases t (no decreasing as it is not large enough!) but the interval is same as the positive interval of randomizing t

                        //checking SHART= just k trains are in the line
                        if (b>k && s==1){
                            Shart= t[b][s]>tt[b-k][S];
                            }
                        }

                    Landa(s,b,t[b][s]);     //computes landa of s and b
                    n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
                    cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
                    if (cond>=n[b][s]){
                        nn[b][s]=0;
                        Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                        }
                    else{
                        nn[b][s]=n[b][s]-cond;
                        Nin[b][s+1]=C;
                        }
                    tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
                    }
            	}
        	}

        //The Non Square Block
        for (int b=b0+k; b<B+1;b++){
            if (b==b0){
                for (int s=s0+1; s<S+1;s++){

                    t[b][s]=max(tt[b-1][s]+tb,tt[b][s-1]+ts);
					//printf("line=388 b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
					
                    //checking SHART= just k trains are in the line
                    if (b>k && s==1){
                        Shart= t[b][s]>tt[b-k][S];
                        //cout<<Shart<<endl;
                        }
                    //check whether the constraints are satisfied of not
                    while (Shart==0){
                        t[b][s]+= q;    //it just increases t 
                        
                        //checking SHART= just k trains are in the line
                        if (b>k && s==1){
                            Shart= t[b][s]>tt[b-k][S];
                            }
                        }

                    Landa(s,b,t[b][s]);     //computes landa of s and b
                    n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
                    cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
                    if (cond>=n[b][s]){
                        nn[b][s]=0;
                        Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                        }
                    else{
                        nn[b][s]=n[b][s]-cond;
                        Nin[b][s+1]=C;
                        }
                    tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
                    }
                }
            else{
                for (int s=1; s<S+1;s++){

                    t[b][s]=max(tt[b-1][s]+tb,tt[b][s-1]+ts);
					//printf("line=423 b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
					
                    //checking SHART= just k trains are in the line
                    if (b>k && s==1){
                        Shart= t[b][s]>tt[b-k][S];
                        
                        }
                    //check whether the constraints are satisfied of not
                    while (Shart==0){
                        t[b][s]+= q;      //it just increases t 
                        
                        //checking SHART= just k trains are in the line
                        if (b>k && s==1){
                            Shart= t[b][s]>tt[b-k][S];
                            }
						}
                        
					
                    Landa(s,b,t[b][s]);     //computes landa of s and b
                    n[b][s]=nn[b-1][s]+landa[b][s]*(t[b][s]-tt[b-1][s]);   //people on station when b arrives s
                    cond=C-Nin[b][s]*(1-beta[s]);   //capacity of the bus
                    if (cond>=n[b][s]){
                        nn[b][s]=0;
                        Nin[b][s+1]=Nin[b][s]*(1-beta[s])+n[b][s];
                        }
                    else{
                        nn[b][s]=n[b][s]-cond;
                        Nin[b][s+1]=C;
                        }
                    tt[b][s]=t[b][s]+tout*beta[s]*Nin[b][s]+tin*(Nin[b][s+1]-(1-beta[s])*Nin[b][s]);  //time of departure
                    }
                }
            }
        }
}



/*The damn main, at last...*/
int main(){
    time(&start);	//time to know how much time each run (for sigle k) takes
    int k=1;	//number of buses in the line
    double alaki=0;		//metropolis dummy variable
	
/*defining files, opening them and writing their header*/
	//The initial schedule
    FILE *init;    //t: time of arrivals
	FILE *initt;	//tt: time of departure
    FILE *inin;    //n: people on the station
    FILE *ininn;   //nn: people who couldn't get in and waiting for another train
    FILE *iniNin;    //N: passengers in the train
    FILE *iniEn;   //Energy
    FILE *inil;    //Landa: the rate of coming passengers to the station (a function of time and station)
    FILE *inib;    //beta: the probability of getting off the train

	init  = fopen ("initial t.txt","w");
    initt = fopen ("initial tt.txt","w");
    inin  = fopen ("initial n.txt","w");
    ininn = fopen ("initial nn.txt","w");
    iniNin  = fopen ("initial Nin.txt","w");
	iniEn = fopen ("initial En.txt","w");
	inil  = fopen ("initial Landa.txt","w");
	inib  = fopen ("initial Beta.txt","w");
	
	
    //The final ones
    FILE *fint;    //t: time of arrivals
	FILE *fintt;	//tt: time of departure
    FILE *finn;    //n: people on the station
    FILE *finnn;   //nn: people who couldn't get in and waiting for another train
    FILE *finNin;    //N: passengers in the train
    FILE *finEn;   //Energy
    FILE *finl;    //Landa: the rate of coming passengers to the station (a function of time and station)
    FILE *finb;    //beta: the probability of getting off the train
	
    fint  = fopen ("final t.txt","w");
    fintt = fopen ("final tt.txt","w");
    finn  = fopen ("final n.txt","w");
    finnn = fopen ("final nn.txt","w");
    finNin  = fopen ("final Nin.txt","w");
	finEn = fopen ("final En.txt","w");
	finl  = fopen ("final Landa.txt","w");
	finb  = fopen ("final Beta.txt","w");
	
    
    //The Important output
    FILE *fileEk;
    fileEk = fopen ("Ek.txt","w");


/*initial values*/
    srand(time(NULL));
    Beta();     //executing beta function


/*k while*/
    while (k<=B){
        //set the system at its initial temperature, reseting all matrix
        memcpy(Nin,zeroarray,sizeof(zeroarray)); //Number of passengers in bus
        memcpy(n,zeroarray,sizeof(zeroarray))  ; //number of passengers in station: in the moment of bus's arrival
        memcpy(nn,zeroarray,sizeof(zeroarray)) ; //number of passengers in station: after leaving the bus
        memcpy(t,zeroarray,sizeof(zeroarray))  ; //schedule
        memcpy(tt,zeroarray,sizeof(zeroarray)) ; //departure schedule
        memcpy(En,zeroarray,sizeof(zeroarray)) ; //Energy matrix


        kT=kTmax;
        E=0;
        Tdet(k);       //making a pre-determined t (with specific t-line headway)
        Energy(1,1,k,tbar);  //evaluating E
        
        
		//Writing initial files
		for (int b=1;b<B+1;b++){
            for (int s=1;s<S+1;s++){
                if (s!=S){
                		fprintf(init,"%f\t",t[b][s]);
                		fprintf(initt,"%f\t",tt[b][s]);
                		fprintf(inin,"%f\t",n[b][s]);
                		fprintf(ininn,"%f\t",nn[b][s]);
                		fprintf(iniNin,"%f\t",Nin[b][s]);
                		fprintf(inil,"%f\t",landa[b][s]);
                		fprintf(inib,"%f\t",beta[s]);
                		fprintf(iniEn,"%f\t",En[b][s]);
                		}
                else if (s==S && b==B){
                		fprintf(init,"%f\n\n\n",t[b][s]);
                		fprintf(initt,"%f\n\n\n",tt[b][s]);
                		fprintf(inin,"%f\n\n\n",n[b][s]);
                		fprintf(ininn,"%f\n\n\n",nn[b][s]);
                		fprintf(iniNin,"%f\n\n\n",Nin[b][s]);
                		fprintf(inil,"%f\n\n\n",landa[b][s]);
                		fprintf(inib,"%f\n\n\n",beta[s]);
                		fprintf(iniEn,"%f\n\n\n",En[b][s]);
                		}
                else{
                		fprintf(init,"%f\n",t[b][s]);
                		fprintf(initt,"%f\n",tt[b][s]);
                		fprintf(inin,"%f\n",n[b][s]);
                		fprintf(ininn,"%f\n",nn[b][s]);
                		fprintf(iniNin,"%f\n",Nin[b][s]);
                		fprintf(inil,"%f\n",landa[b][s]);
                		fprintf(inib,"%f\n",beta[s]);
                		fprintf(iniEn,"%f\n",En[b][s]);
                	}
                }//end of s
            }// end of b
            
		
/*kT while*/
		while (kT>1.5){
			
            r=0;	//resetting acceptance ratio
            
/*main loops for b and s, each imax times*/
            for (int b=1;b<B+1;b++){	//for all buses
                for (int s=1;s<S+1;s++){	//for all stations
                    for (int i=0;i<imax;i++){	//randomize t for imax times
						
                        //copying into initials
                        memcpy(Nini,Nin,sizeof(Nin));
                        memcpy(ti,t,sizeof(ti));
                        memcpy(tti,tt,sizeof(tti));
                        memcpy(ni,n,sizeof(ni));
                        memcpy(nni,nn,sizeof(nni));
                        memcpy(Eni,En,sizeof(Eni));
                        memcpy(landai,landa,sizeof(landai));
                        
						//randomizing and evaluating energy
						Ei=En[B][S];	//initial energy
                        Func(b,s,k);

                        Energy(b,s,k,tline);
                        Ef=En[B][S];	//final energy
                        dE=Ef-Ei;		//dE
                        
                        //saving Emin
                        if (Ef-Emin<0){
                            Emin=Ef;
                        	}

                        alaki=(double)rand()/(double)RAND_MAX;	//metropolis random parameter
						
						//metripolis condition
                        if (alaki>=exp(double(-dE/kT))){	//when it is not accepted!
                        	//roll back all to their initials
                        	
                            memcpy(Nin,Nini,sizeof(Nin));
                            memcpy(t,ti,sizeof(t));
                            memcpy(tt,tti,sizeof(tt));
                            memcpy(n,ni,sizeof(n));
                            memcpy(nn,nni,sizeof(nn));
                            memcpy(En,Eni,sizeof(En));
                            memcpy(landa,landai,sizeof(landa));
                            }
                        else{	//when it is accepted!
                        	
                            //cout<<"dE= "<<dE<<"\t E="<<Ef<<"\t kT="<<kT<< endl;
                            r++;
                            //cout<<t[3][1]<<"\tkhbadlkfhadjkhakf\t"<<tt[2][S]<<endl;
							}
						
                		}// i for end
            		}//s for end
        		}//b for end
        		
            R=r/double(imax*B*S);	//normalizing r to R
            
            //output foe myself
			//cout<<"r="<<R<<"\tkT="<<kT<<"\tk="<<k<<"\tE="<<En[B][S]<<"\tdelta="<<delta<<endl;

			//controling delta
            if (R<0.5){
                delta*=0.5;
                }
            else{
                delta*=2;
                if (delta>T/2){
                	delta=T/2;
					}
                }
            
            kT*=0.9;	//colling the system
			}//end of mail loop while for kT

        //time(&sample);
        //cout<<"\ttime (seconds)="<<long(sample-start)<<endl;
		fprintf(fileEk,"%d\t%f\t%f\n",k,En[B][S],Emin);	//write E of k and minimum energy of system
        
		
        //Writing final files
		for (int b=1;b<B+1;b++){
            for (int s=1;s<S+1;s++){
                if (s!=S){
                	//printf("b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
                	//printf("b= %d and s=%d and k=%d and tt[b][s]=%f\n",b,s,k,tt[b][s]);
					fprintf(fint,"%f\t",t[b][s]);
                	fprintf(fintt,"%f\t",tt[b][s]);
                	fprintf(finn,"%f\t",n[b][s]);
                	fprintf(finnn,"%f\t",nn[b][s]);
                	fprintf(finNin,"%f\t",Nin[b][s]);
                	fprintf(finl,"%f\t",landa[b][s]);
                	fprintf(finb,"%f\t",beta[s]);
                	fprintf(finEn,"%f\t",En[b][s]);
                	}
                else{
				 if (s==S && b==B){
				 	//printf("b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
                	//printf("b= %d and s=%d and k=%d and tt[b][s]=%f\n",b,s,k,tt[b][s]);
                	fprintf(fint,"%f\n\n\n",t[b][s]);
                	fprintf(fintt,"%f\n\n\n",tt[b][s]);
                	fprintf(finn,"%f\n\n\n",n[b][s]);
                	fprintf(finnn,"%f\n\n\n",nn[b][s]);
                	fprintf(finNin,"%f\n\n\n",Nin[b][s]);
                	fprintf(finl,"%f\n\n\n",landa[b][s]);
                	fprintf(finb,"%f\n\n\n",beta[s]);
                	fprintf(finEn,"%f\n\n\n",En[b][s]);
                	}
                else{
                	
                	//cout<<tt[b][S]<<"\t"<<t[b+1][1]<<endl;
                	//printf("b= %d and s=%d and k=%d and t[b][s]=%f\n",b,s,k,t[b][s]);
                	//printf("b= %d and s=%d and k=%d and tt[b][s]=%f\n",b,s,k,tt[b][s]);
                	fprintf(fint,"%f\n",t[b][s]);
                	fprintf(fintt,"%f\n",tt[b][s]);
                	fprintf(finn,"%f\n",n[b][s]);
                	fprintf(finnn,"%f\n",nn[b][s]);
                	fprintf(finNin,"%f\n",Nin[b][s]);
                	fprintf(finl,"%f\n",landa[b][s]);
                	fprintf(finb,"%f\n",beta[s]);
                	fprintf(finEn,"%f\n",En[b][s]);
                	}
                	}
                }//end of s
            }// end of b
            
        k++;	//add one bus to the system and run again

    	}

cout << '\a';	//ring when it is finished
}
