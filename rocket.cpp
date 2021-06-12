#include <vector>
#include <tuple>
#include <iostream>
#include "rocket.h"
#include <cmath>
#include <iterator>

typedef std::vector<double> (*vec2vec)(std::vector<double>);
//DEFINING FUNCTIONS OF THE BODY CLASS

//Constructor
body::body(std::vector<double> initState,double mass,const int dim,const int id){
    this->state = initState;
    this->M = mass;
    this->dim = dim;
    this->trajectory.push_back(initState);
    this->id = id;
}

//Updates the internal state
void body::update(bool save=false){
    this->state = this->nextstate;
    if(save ==true){
        trajectory.push_back(this->nextstate);  //Probably we only need to store a smaller number of coordinates than is generated
    }
}
void body::update(double newMass){
    this->M = newMass;
}

//Some getter functions
std::vector<double> body::getR(){
    return std::vector<double> {this->state[1],this->state[2]};}

std::vector<double> body::getV(){
    return std::vector<double> {this->state[3],this->state[4]};
}
double body::Ek(){
    double vSq = 0;
    for(int i=0; i<=sizeof(this->getV()); i++){
        vSq+= this->getR()[i]*this->getV()[i];
    }
    return 0.5*this->M*vSq;
}

//Updates force on body due to a second body b
void body::force(std::vector<double> *state,body *b,int dim,std::vector<double> *force){
    //Computes magnitude of the gravitational force on a body due to a body b
    //Evaluates this and returns a double
    double rSq = 0;
    for(int i=0; i<dim;i++){
        rSq += pow((*state)[i]-b->getR()[i],2);
    }
    for(int i=0; i<dim;i++){      
        double val;
        //(*force)[i]=1;
        val = -1*(b->M)*((*state)[i]-b->getR()[i])/pow(rSq,1.5);
        (*force)[i] = val;
        //(*force)[i]+=-b->M*((*state)[1+i]-b->getR()[1+i])/pow(rSq,1.5);
    }
}

void body::RK4(std::vector<body*> bodies,double dt){
    //Uses RK4 method to generate statepoint at next timestep
    std::vector<double> istate;               //Intermediate RK4 state
    for(int k=1;k<this->state.size();k++){
        istate.push_back(this->state[k]);
    }
    
    std::vector<double>  fstate;              //The updated final state
    std::vector<double> *iptr = & istate;     //Points to intermediate RK4 states
    std::vector<double*> fCalls;              //Vector of pointers to the result of function calls used to update fstate
    int step = 1;
    
    std::vector<double> cforce{0,0};        //To store the cumulative force on the body
    std::vector<double> *f;
    f = &cforce;
    
    for(int i =0; i<2*dim;i++){             //Pointers to the numbers which are used to update the state
        double *p;
        if(i<dim){
            p = & istate[dim+i];}
        else{
            p = & cforce[i-dim];}
        fCalls.push_back(p);
    }
    
    while(step<=4){
        *f = {0,0};
        //Initialize final state to to be the initial state
        if(step==1){
            for(int c=0;c<this->state.size();c++){
            fstate.push_back((this->state)[c]);}
        }

        //Calculate cumulative force on object
        for(int i=0;i<bodies.size();i++){
            if(i!=this->id){                    //Prevents unphysical self-interaction
            force(iptr,bodies[i],this->dim,f);} //Updates cforce by reference
        }

        //Update final state using function calls
        for(int c=0;c<istate.size();c++){
            if(step==1||step==4){
                fstate[c+1] += *(fCalls[c])*dt/6;}
            else{
                fstate[c+1] += *(fCalls[c])*dt/3;}
        }

        //Update internal state
        if(step<3){
            for(int c=0;c<istate.size();c++){
                if(c<dim){
                    istate[c]+=dt*istate[c+dim]/2;}
                else{
                    istate[c] +=dt*cforce[c-dim]/2;}   
            }
        }
        else{
            for(int c=0;c<istate.size();c++){
                if(c<dim){
                    istate[c]+=dt*istate[c+dim];} 
                else{
                    istate[c] +=dt*cforce[c-dim];}    
            }
        }
    step++;
    }
    fstate[0]+=dt;
    this->nextstate = fstate;
}
int main(){
    int nIters = 1000;
    int dim = 2;

    double mass = 1.0;
    double timestep = pow(10,-2);
    auto s1 = std::vector<double> {0,0,0,0,0};       //Fixed mass
    auto s2 = std::vector<double> {0,0,7,0.36,0};    //Orbiter. With same mass objects these parameters give a fairly stable elliptic orbit
    
    //body earth(s2,mass);
    std::vector<body*> planets;
    planets.push_back(new body(s1,mass,dim,0));
    planets.push_back(new body(s2,mass,dim,1));
    
    for(int i=0;i<1;i++){
    planets[1]->RK4(planets,timestep);
    planets[1]->update();
    }
    
    for(int k=0;k<5;k++){       //Printing out the state after one timestep to see if it agrees with the Python implementation
        std::cout<<(planets[1]->state)[k]<<std::endl;
    }   
    /*
    THIS CODE WOULD BE USED FOR THE GENERAL N-BODY SIMULATION
    for(int i=0;i<nIters;i++){
        for(int k = 0;k<planets.size();k++){
            planets[k]->RK4(planets,timestep);  //RK4 update to find next states
        }
        for(int k = 0;k<planets.size();k++){
            if(i%10==0){
            planets[k]->update(true);}         //Save trajectory each 10 iterations
            else{
            planets[k]->update();
            }
        }
    }
    */
    return 0;
}
