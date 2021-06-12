#ifndef _body
#define _body

//Massive body class
typedef std::vector<double> (*vec2vec)(std::vector<double>);
class body{
    public:
        //DATA
        std::vector<double> state;                  //State of system
        double M;                                   //Mass
        //FUNCTIONS
        body(std::vector<double> initState,double mass,int dim,int id); //Constructor
        std::vector<double> getR();                 //Return position
        std::vector<double> getV();                 //Return momentum
        double Ek();                                //Calculate kinetic energy
        void force(std::vector<double> *pos,body *b,int dim,std::vector<double> *force);    //Calculate force on body due to another
        void RK4(std::vector<body*> bodies,double dt);  //Performs one RK4 update
        void update(bool save);    //Updates state using the result of RK4 which is stored in a private variable until update
        void update(double mass);  //Overloaded defintion of update which instead allows for variable mass 
    private:
        int id;     //Each object should have a unique identifier
        int dim;    //Dimentionality of problem; dim=2 for in plane motion                              
        std::vector<double> nextstate;                  //Result of RK4 stored here until update called
        std::vector<std::vector<double>> trajectory;   //Trajectory of the particle
};
#endif
