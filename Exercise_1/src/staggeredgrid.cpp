#include "staggeredgrid.h"
#include <iostream>

StaggeredGrid::StaggeredGrid(std::array<int,2> size) 
:pressure(size), velocity_X(size), velocity_Y(size), F(size), G(size), RHS(size)
{
}

StaggeredGrid::~StaggeredGrid()
{
}
// hab vergessen, was wir an den Ecken nehmen, muss dann noch angepasst werden
void StaggeredGrid::setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2> right,std::array<double,2> bottom)
{  int i_max = pressure.size()[0], j_max = pressure.size()[1];
   std::cout<<"Wertzuweisung an den Ecken anpassen"<<std::endl; 
    for (int j = 0; j < j_max; j++)
   {
        velocity_X(0,j)=top[0];
        velocity_Y(0,j)=top[1];   
    
        velocity_X(j_max-1,j)=bottom[0];
        velocity_Y(j_max-1,j)=bottom[1]; 
   }
   
   for (int i = 0; i < i_max; i++)
   {
        velocity_X(i,0)=left[0];
        velocity_Y(i,0)=left[1];  

        velocity_X(i,i_max-1)=right[0];
        velocity_Y(i,i_max-1)=right[1]; 
        
   }
      
}

void StaggeredGrid::print(std::string str)
{
    if (str == "pressure")
    {
        pressure.print();
    } else if (str == "velocity_X" || str == "velocity_x")
    {
        velocity_X.print();
    }else if (str == "velocity_Y" || str == "velocity_y")
    {
        velocity_Y.print();
    }else if (str == "F")
    {
        F.print();
    }else if (str == "G")
    {
        G.print();
    }else if (str == "all")
    {
        std::cout<< "pressure" << std::endl;   
        pressure.print();
        std::cout<< "velocity_X" << std::endl;  
        velocity_X.print();
        std::cout<< "velocity_Y" << std::endl;  
        velocity_Y.print();
        std::cout<< "F" << std::endl;  
        F.print();
        std::cout<< "F" << std::endl;  
        G.print();
    }





}

std::array<float,2> StaggeredGrid::meshWidth()  const
{
    //TODO
}
std::array<int,2> StaggeredGrid::nCells() const
{
    //TODO
}
FieldVariable StaggeredGrid::p() const
{
    //TODO
}
float StaggeredGrid::p(int i, int j) const
{
    //TODO
}
FieldVariable StaggeredGrid::u() const
{
    //TODO
}
float StaggeredGrid::u(int i, int j) const
{
    //TODO
}
FieldVariable StaggeredGrid::v() const
{
    //TODO
}
float StaggeredGrid::v(int i, int j) const
{
    //TODO
}
float StaggeredGrid::dx() const
{
    //TODO
}
float StaggeredGrid::dy() const
{
    //TODO
}
int StaggeredGrid::uIBegin() const
{
    //TODO
} 
int StaggeredGrid::uIEnd() const
{
    //TODO
} 
int StaggeredGrid::uJBegin() const
{
    //TODO
} 
int StaggeredGrid::uJEnd() const
{
    //TODO
} 
int StaggeredGrid::vIBegin() const
{
    //TODO
} 
int StaggeredGrid::vIEnd() const
{
    //TODO
} 
int StaggeredGrid::vJBegin() const
{
    //TODO
} 
int StaggeredGrid::vJEnd() const
{
    //TODO
} 
int StaggeredGrid::pIBegin() const
{
    //TODO
} 
int StaggeredGrid::pIEnd() const
{
    //TODO
} 
int StaggeredGrid::pJBegin() const
{
    //TODO
} 
int StaggeredGrid::pJEnd() const
{
    //TODO
} 