#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"

//constructor
Discretization::Discretization(Settings settings):
 StaggeredGrid(settings.nCells) //, theGrid(settings.nCells)
{
setBorderVelocity(settings.dirichletBcTop,settings.dirichletBcLeft,settings.dirichletBcRight,settings.dirichletBcBottom);
settings_=settings;
}

//destructor
Discretization::~Discretization()
{
}


void Discretization::calculation()
{
    
std::cout<< "u und v durch direkten Zugriff auf die FieldVariablen ersetzen" <<std::endl;
std::cout<< "wenn debugging erledigt ist UND DT anpassen" <<std::endl;
FieldVariable u=velocity_X;
FieldVariable v=velocity_Y;
float delta_x                 = settings_.physicalSize[0] / settings_.nCells[0];
float delta_y                 = settings_.physicalSize[1] / settings_.nCells[1];
float deltat                  = min3(min2(delta_x,delta_y)/4,delta_x/v.absmax(),delta_y/v.absmax())*settings_.tau;
float delta_x_divisor         = 1 / (delta_x);
float delta_y_divisor         = 1 / (delta_y);
float delta_x_quadrat_divisor = 1 / (delta_x * delta_x);
float delta_y_quadrat_divisor = 1 / (delta_y * delta_y);

//vordefinieren der speicherskalare
float u_xx                  ;
float u_yy                  ;
float v_xx                  ;
float v_yy                  ;
float u_quadrat_x           ;
float v_quadrat_y           ;
float uv_x                  ;
float uv_y                  ;
float u_x                   ;
float v_y                   ;

float zwischenwert_uiphj    ; 
float zwischenwert_uimhj    ;
float zwischenwert_vijph    ;
float zwischenwert_vijmh    ;
float zwischenwert_viphj    ;
float zwischenwert_viphjm1  ;
float zwischenwert_uijph    ;    
float zwischenwert_uijmh    ;
float zwischenwert_vimhj    ;
float zwischenwert_uim1jph  ;

for (int j = 1; j < settings_.nCells[1]-1; j++)
{
    for (int i = 1; i < settings_.nCells[0]-1; i++)
    {
        zwischenwert_viphj    = v(i+1,j)   + v(i,j)   ;
        zwischenwert_viphjm1  = v(i+1,j-1) + v(i,j-1) ;
        zwischenwert_uijph    = u(i,j+1)   + u(i,j)   ; // kann ersetzt werden
        zwischenwert_uijmh    = u(i,j)     + u(i,j-1) ;
        zwischenwert_vimhj    = v(i,j)     + v(i-1,j) ;
        zwischenwert_uim1jph  = u(i-1,j+1) + u(i-1,j) ;

        zwischenwert_uiphj = u(i+1,j) + u(i,j)   * u(i+1,j) + u(i,j)   ;
        zwischenwert_uimhj = u(i,j)   + u(i-1,j) * u(i,j)   + u(i-1,j) ;
        zwischenwert_vijph = u(i,j+1) + u(i,j)   * u(i,j+1) + u(i,j)   ;
        zwischenwert_vijmh = u(i,j)   + u(i,j-1) * u(i,j)   + u(i,j-1) ;

        u_xx=delta_x_quadrat_divisor * (u(i+1,j) - 2 * u(i,j) + u(i-1,j));
        u_yy=delta_y_quadrat_divisor * (u(i,j+1) - 2 * u(i,j) + u(i,j-1));

        v_xx = delta_x_quadrat_divisor * (v(i+1,j) - 2 * v(i,j) + v(i-1,j));
        v_yy = delta_y_quadrat_divisor * (v(i,j+1) - 2 * v(i,j) + v(i,j-1));

        u_quadrat_x = delta_x_divisor / 4 * (zwischenwert_uiphj - zwischenwert_uimhj);
        v_quadrat_y = delta_y_divisor / 4 * (zwischenwert_vijph - zwischenwert_vijmh);

        uv_x = delta_x_divisor /4 * (zwischenwert_viphj * zwischenwert_uijph - zwischenwert_vimhj   * zwischenwert_uim1jph);
        uv_y = delta_y_divisor /4 * (zwischenwert_viphj * zwischenwert_uijph - zwischenwert_viphjm1 * zwischenwert_uijmh);

        u_x = delta_x_divisor * (u(i,j)-u(i-1,j));
        v_y = delta_y_divisor * (v(i,j)-v(i,j-1));

        F(i,j) = u(i,j) + deltat * ((u_xx + u_yy) / settings_.re - uv_y - u_quadrat_x + settings_.g[0]);
        G(i,j) = v(i,j) + deltat * ((v_xx + v_yy) / settings_.re - uv_x - v_quadrat_y + settings_.g[1]);
    }
}
}


float Discretization::min2(float value1, float value2) const
{
    if (value1<=value2)
    {
        return value1;
    }else 
    {
        return value2;
    }
}

float Discretization::min3(float value1, float value2, float value3) const
{
    if (value1<=value2 && value1<=value3)
    {
        return value1;
    }else if (value2<=value1 && value2<=value3)
    {
        return value2;
    } else
    {
        return value3;
    }
}

float Discretization::f(int i, int j) const
{
    //TODO
}
float Discretization::g(int i, int j) const
{
    //TODO
}
float Discretization::rhs(int i, int j) const
{
    //TODO
}
