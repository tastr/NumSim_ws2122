#pragma once

#include "fieldvariable.h"

class StaggeredGrid 
{
private:
    
public:
    FieldVariable pressure, velocity_X, velocity_Y,F, G; //Sollten eigentlich private sein?!
    StaggeredGrid(std::array<int,2> size);
    ~StaggeredGrid();
     void setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2>  right,std::array<double,2> bottom); 
     void print(std::string str);
     // set and get functions to acess the fieldvariables

     // get functions:
     std::array<float,2> meshWidth(); //TODO
     std::array<int,2> nCells(); //TODO
     FieldVariable p(); //TODO
     float p(int i, int j); //TODO
     FieldVariable u(); //TODO
     float u(int i, int j); //TODO
     FieldVariable v(); //TODO
     float v(int i, int j); //TODO
     float dx(); //TODO
     float dy(); //TODO
     int uIBegin(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int uIEnd(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int uJBegin(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int uJEnd(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vIBegin(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vIEnd(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vJBegin(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vJEnd(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pIBegin(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pIEnd(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pJBegin(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pJEnd(); //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     
    

     //void getValue(std::string str,int i, int j);
};


