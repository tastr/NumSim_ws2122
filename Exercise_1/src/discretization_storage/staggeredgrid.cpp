#include "staggeredgrid.h"
#include "settings.h"
#include <iostream>
#include <cassert>

//StaggeredGrid::StaggeredGrid(std::array<int,2> size) 
StaggeredGrid::StaggeredGrid(Settings settings) 
:pressure({settings.nCells[0]+2,settings.nCells[1]+2})
,velocity_X({settings.nCells[0]+1,settings.nCells[1]+2})
,velocity_Y({settings.nCells[0]+2,settings.nCells[1]+1})
,settings_(settings)
{ 
    setSize_(settings.nCells);
    delta_x=settings_.physicalSize[0] / settings_.nCells[0];
    delta_y=settings_.physicalSize[1] / settings_.nCells[1];
    epsilon=settings.epsilon;
}


void StaggeredGrid::setSize_(std::array<int,2> nCells)
{
size_={nCells[0]+2,nCells[1]+2};
}

std::array<int,2> StaggeredGrid::getSize() const
{
return size_;
}

// Sets the Dirichlet border conditions for the velocity variables
void StaggeredGrid::setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2> right,std::array<double,2> bottom)
{  
//     int i_max = pressure.size()[0], j_max = pressure.size()[1]; 
//     for (int j = 0; j < j_max; j++)
//    {
//         velocity_X(0,j)=left[0];
//         velocity_Y(0,j)=2*left[1]-velocity_Y(1,j);   
    
//         velocity_X(j_max-1,j)=right[0];
//         velocity_Y(j_max-1,j)=2*right[1]-velocity_Y(j_max-2,j); 
//    }
   
//    // i starts at 1 and goes to i_max-1 so that the wall is the BC in corners
//    for (int i = 1; i < i_max-1; i++)
//    {
//         velocity_X(i,0)=2*bottom[0]-velocity_X(i,1);
//         velocity_Y(i,0)=bottom[1];  

//         velocity_X(i,i_max-1)=2*top[0]-velocity_X(i,i_max-2);
//         velocity_Y(i,i_max-1)=top[1]; 
        
//    }

   // set u velocity
   int i_u_max = velocity_X.size()[0], j_u_max = velocity_X.size()[1];
   for (int j = 0; j < j_u_max; j++)
   {
       velocity_X(0,j)=left[0];
        velocity_X(i_u_max-1,j)=right[0];
   }
   for (int i = 1; i < i_u_max-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
   {
       velocity_X(i,0)=2*bottom[0]-velocity_X(i,1);
       velocity_X(i,j_u_max-1)=2*top[0]-velocity_X(i,j_u_max-2);
   }

   // set v velocity
   for (int j = vJBegin(); j < vJEnd(); j++)
   {
       velocity_Y(0,j)=2*left[1]-velocity_Y(1,j); 
       velocity_Y(vIEnd()-1,j)=2*right[1]-velocity_Y(vIEnd()-2,j);

   }
   for (int i = vIBegin()+1; i < vIEnd()-1; i++)
   {
       velocity_Y(i,0)=bottom[1]; 
       velocity_Y(i,vJEnd()-1)=top[1]; 
   }
   
   
   
    

   
      
}

//updates the border terms for the pressure variable in the end of each iteration
void StaggeredGrid::updatedPressureBC()
{
    int i_max = pressure.size()[0], j_max = pressure.size()[1];
    for (int j = 0; j < j_max; j++)
    {
        pressure(0,j)=pressure(1,j);
        pressure(i_max-1,j)=pressure(i_max-2,j);
    }
    for (int i = 0; i < i_max; i++)
    {
        pressure(i,0)=pressure(i,1);
        pressure(i,j_max-1)=pressure(i,j_max-2);
    }
}

//prints out the values stored in the FieldVariables
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
    }else if (str == "all")
      {
        std::cout<< "pressure" << std::endl;   
        pressure.print();
        std::cout<< "velocity_X" << std::endl;  
        velocity_X.print();
        std::cout<< "velocity_Y" << std::endl;  
        velocity_Y.print();
        
      }





}

// functions return values describing the grid
std::array<double,2> StaggeredGrid::meshWidth()  const
{
 return {dx(),dy()}  ; 
}
std::array<int,2> StaggeredGrid::nCells() const
{
    return settings_.nCells;
}
double StaggeredGrid::getepsilon() const
{
 return epsilon;
}
double StaggeredGrid::dx() const
{
    return delta_x;
}
double StaggeredGrid::dy() const
{
    return delta_y;
}


//functions to acess entire Fieldvariables of the grid
FieldVariable StaggeredGrid::p() const
{
    return pressure;
}
double StaggeredGrid::p(int i, int j) const
{
    return pressure(i,j);
}
FieldVariable StaggeredGrid::u() const
{
    return velocity_X;
}


//functions to retrieve singular values from the FieldVariables
double StaggeredGrid::u(int i, int j) const
{
    return velocity_X(i,j);
}
FieldVariable StaggeredGrid::v() const
{
    return velocity_Y;
}
double StaggeredGrid::v(int i, int j) const
{
    return velocity_Y(i,j);
}


//functions return the size values of the FieldVariables 
int StaggeredGrid::uIBegin() const
{
    return 0;
} 
int StaggeredGrid::uIEnd() const
{
    return velocity_X.size()[0];
} 
int StaggeredGrid::uJBegin() const
{
    return 0;
} 
int StaggeredGrid::uJEnd() const
{
    return velocity_X.size()[1];
} 
int StaggeredGrid::vIBegin() const
{
    return 0;
} 
int StaggeredGrid::vIEnd() const
{
    return velocity_Y.size()[0];
} 
int StaggeredGrid::vJBegin() const
{
    return 0;
} 
int StaggeredGrid::vJEnd() const
{
    return velocity_Y.size()[1];
} 
int StaggeredGrid::pIBegin() const
{
    return 0;
} 
int StaggeredGrid::pIEnd() const
{
    return pressure.size()[0];
} 
int StaggeredGrid::pJBegin() const
{
    return 0;
} 
int StaggeredGrid::pJEnd() const
{
    return pressure.size()[1];
} 

//functions to write values into the Fieldvariables 
void StaggeredGrid::setU(int i, int j,double value)
{
velocity_X(i,j)=value;
}
void StaggeredGrid::setV(int i, int j,double value)
{
velocity_Y(i,j)=value;
}
void StaggeredGrid::setP(int i, int j,double value)
{
pressure(i,j)=value;
}


// functions to overwrite the FieldVariables of the grid with the Values of other Fieldvariables with maching size
void StaggeredGrid::setV(FieldVariable value)
{  // assert that indices are in range
assert(value.size()[0] == velocity_Y.size()[0]);
assert(value.size()[1] == velocity_Y.size()[1]);
for (int j = 0; j < value.size()[1] ; j++)
{
    for (int i = 0; i < value.size()[0] ; i++)
    {
        velocity_Y(i,j)=value(i,j);  
    }
}
}

void StaggeredGrid::setP(FieldVariable value)
    {  // assert that indices are in range
    assert(value.size()[0] == pressure.size()[0]);
    assert(value.size()[1] == pressure.size()[1]);
    for (int j = 0; j < value.size()[1] ; j++)
    {
        for (int i = 0; i < value.size()[0] ; i++)
        {
            pressure(i,j)=value(i,j);  
        }
    }
}


//function to calculate the absolute value of a double
double StaggeredGrid::abs(double number) const
   {
 if (number>=0)
  {
    return number;
  }else
  {
    return -number;
  }
   }


