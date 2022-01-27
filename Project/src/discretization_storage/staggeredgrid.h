#pragma once

#include "fieldvariable.h"
#include "settings.h"
#include <cassert>

class StaggeredGrid 
{
protected:
    Settings settings_;
    std::array<int,2> size_;
    FieldVariable pressure,velocity_X , velocity_Y, type;
    void setSize_(std::array<int,2> nCells);
    double delta_x;
    double delta_y;
    double epsilon;
public:
    
    //StaggeredGrid(std::array<int,2> size);    ~StaggeredGrid();
     StaggeredGrid(Settings settings);
     void setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2>  right,std::array<double,2> bottom); 
     void updatedPressureBC();
     void print(std::string str);
     double abs(double number) const;
     
     // set and get functions to acess the fieldvariables

     // get functions:
     // werden gebraucht für den output_writer
     std::array<double,2> meshWidth() const; //Hab mal  xmax-xmin bzw ymax-ymin eingefuegt
     std::array<int,2> getSize() const; 
     
     std::array<int,2> nCells() const; 
     FieldVariable p() const; 
     double p(int i, int j) const; 
     FieldVariable u() const; 
     double u(int i, int j) const; 
     FieldVariable v() const; 
     double v(int i, int j) const;
     //FieldVariable typ() const; 
     double typ(int i, int j) const;

     double dx() const; //ist das nicht identisch zu meshwidth?
     double dy() const; //ist das nicht identisch zu meshwidth?
     double getepsilon() const; //ist das nicht identisch zu meshwidth?

     int getMaxIteration() const;
     
     int uIBegin() const;  
     int uIEnd() const;  
     int uJBegin() const;  
     int uJEnd() const;  
     int vIBegin() const;  
     int vIEnd() const;  
     int vJBegin() const;  
     int vJEnd() const;  
     int pIBegin() const;  
     int pIEnd() const;  
     int pJBegin() const;  
     int pJEnd() const;  
     
     // get function fuer FieldVariablen
     void setU(int i, int j, double);
     void setV(int i, int j, double);
     void setP(int i, int j, double);
         
     // Keine Ahnung ob man die braucht, aber ne function um die FV komplett zu ersetzen 
     void setU(FieldVariable value);
     void setV(FieldVariable value);
     void setP(FieldVariable value);    
         
     void setObstacleVelocity();      
     void setObstaclePressure(int i, int j);      

     //void getValue(std::string str,int i, int j);
};


