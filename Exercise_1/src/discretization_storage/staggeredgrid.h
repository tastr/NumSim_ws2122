#pragma once

#include "fieldvariable.h"
#include "settings.h"
#include <cassert>

class StaggeredGrid 
{
protected:
    Settings settings_;
    std::array<int,2> size_;
    FieldVariable pressure,velocity_X , velocity_Y;
    void setSize_(std::array<int,2> nCells);
    double delta_x;
    double delta_y;
    double epsilon;
public:
    
    //StaggeredGrid(std::array<int,2> size);    ~StaggeredGrid();
     StaggeredGrid(Settings settings);
     void setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2>  right,std::array<double,2> bottom); 
     void print(std::string str);
     double abs(double number) const;
     
     // set and get functions to acess the fieldvariables

     // get functions:
     // werden gebraucht für den output_writer
     std::array<double,2> meshWidth() const; //Hab mal  xmax-xmin bzw ymax-ymin eingefuegt
     std::array<int,2> getSize() const; 
     
     std::array<int,2> nCells() const; //TODO
     FieldVariable p() const; //TODO
     double p(int i, int j) const; //TODO
     FieldVariable u() const; //TODO
     double u(int i, int j) const; //TODO
     FieldVariable v() const; //TODO
     double v(int i, int j) const; //TODO
     double dx() const; //ist das nicht identisch zu meshwidth?
     double dy() const; //ist das nicht identisch zu meshwidth?
     double getepsilon() const; //ist das nicht identisch zu meshwidth?
     
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
     
     // get function fuer FieldVariablen
     void setU(int i, int j, double);
     void setV(int i, int j, double);
     void setP(int i, int j, double);
         
     // Keine Ahnung ob man die braucht, aber ne function um die FV komplett zu ersetzen 
     void setU(FieldVariable value);
     void setV(FieldVariable value);
     void setP(FieldVariable value);    
         
         

     //void getValue(std::string str,int i, int j);
};


