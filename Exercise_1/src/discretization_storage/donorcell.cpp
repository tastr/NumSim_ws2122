#include "donorcell.h"
#include "discretization.h"



   DonorCell::DonorCell(Settings settings):Discretization(settings) 
   {
   }

   DonorCell::~DonorCell() 
   {
   }


(u(i,j) + u(i+1,j)) * (u(i,j) + u(i+1,j))
(u(i-1,j) + u(i,j)) * (u(i-1,j) + u(i,j))
abs(u(i,j)+u(i+1,j)) * (u(i,j)-u(i+1,j))
abs(u(i-1,j)+u(i,j)) * (u(i-1,j)-u(i,j))

uquadratX = (  ) )/(2*dx());


(v(i,j)+v(i+1,j)) * (u(i,j)+u(i,j+1))
(v(i,j-1)+v(i+1,j-1))* (v(i,j-1)+v(i+1,j-1))
abs(v(i,j)+v(i+1,j)) * (u(i,j)-u(i,j+1))
abs(v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)-u(i,j))

uv_y 



(u(i,j) + u(i,j+1)) * (u(i,j) + u(i,j+1))
(u(i,j-1) + u(i,j)) * (u(i,j) + u(i,j-1))
abs(u(i,j)+u(i,j+1)) * (u(i,j)-u(i,j+1))
abs(u(i,j-1)+u(i,j)) * (u(i,j-1)-u(i,j))

uquadratY = (  ) )/(2*dx());


(v(i+1,j)-v(i,j))*abs(u(i,j+1)+u(i,j))
(v(i,j)-v(i-1,j))*abs(u(i-1,j+1) + u(i-1,j)) 



uv_y







