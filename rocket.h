#ifndef _body
#define _body

/**Massive body class
N.B. I have written setter methods instead of a constructor like body::body()
because I want to initialize arrays of objects without setting initial states
which should all be different.Ideally, I would have the functions which
update the internal state as private functions and only call the constructor
once
**/
typedef std::vector<double> (*vec2vec)(std::vector<double>);
class body{
    public:
        //DATA
        std::vector<double> state;                  //State of system
        double M;
        //FUNCTIONS
        body(std::vector<double> initState,double mass,int dim,int id);
        std::vector<double> getR();                 //Return position
        std::vector<double> getV();                 //Return momentum
        double Ek();                                //Calculate kinetic energy
        std::vector<double> propagate(float dt);  //RK4 method for central force problem. Propagates by a timestep dt
        void force(std::vector<double> *pos,body *b,int dim,std::vector<double> *force);//Calculate force on body due to another
        void RK4(std::vector<body*> bodies,double dt);
        void update(bool save);    //Updates state
        void update(double mass);   
    private:
        int id;
        int dim;                                  
        std::vector<double> nextstate;                  //Updates mass
        std::vector<std::vector<double>> trajectory;   //Trajectory of the particle
};

//This lets me add forces together and return an evaluatable function object
std::function<int(int,int)> addForces(int (*const fPtr)(int x),int (*const gPtr)(int y));
void printPos(body bodies);
void printPos_single(body *B);
#endif