#pragma once

#include "fieldvariable.h"

class StaggeredGrid 
{
protected:
    FieldVariable pressure, velocity_X, velocity_Y,F, G; //Sollten eigentlich private sein?!
public:
<<<<<<< HEAD
    FieldVariable pressure, velocity_X, velocity_Y,F, G, RHS;
=======
    
>>>>>>> b58dd6ee07f3c0424f80da1333cd28f86fbe51bf
    StaggeredGrid(std::array<int,2> size);
    ~StaggeredGrid();
     void setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2>  right,std::array<double,2> bottom); 
     void print(std::string str);
     // set and get functions to acess the fieldvariables

     // get functions:
     // werden gebraucht für den output_writer
     std::array<float,2> meshWidth() const; //TODO
     std::array<int,2> nCells() const; //TODO
     FieldVariable p() const; //TODO
     float p(int i, int j) const; //TODO
     FieldVariable u() const; //TODO
     float u(int i, int j) const; //TODO
     FieldVariable v() const; //TODO
     float v(int i, int j) const; //TODO
     float dx() const; //TODO
     float dy() const; //TODO
     int uIBegin() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int uIEnd() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int uJBegin() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int uJEnd() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vIBegin() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vIEnd() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vJBegin() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int vJEnd() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pIBegin() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pIEnd() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pJBegin() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     int pJEnd() const; //TODO keine ahnung was die Funktion machen soll. Vielleicht startindex zurückgeben?
     
    

     //void getValue(std::string str,int i, int j);
     
    
};


