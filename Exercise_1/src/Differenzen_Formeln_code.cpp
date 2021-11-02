//pressure_x(i,j) = delta_x_divisor*(pressure(i+1,j)-pressure(i,j));
//pressure_x(i,j) = delta_y_divisor*(pressure(i,j+1)-pressure(i,j+1));

// input Arrays u und v sowie deltax deltay 

//vordefinieren der oefter verwendeten skalare

float delta_x_divisor         = 1 / delta_x;
float delta_y_divisor         = 1 / delta_y;
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

for (j = 1; j < Tiefe; j++)
{
    for ( i = 1; i < Laenge; i++)
    {
        zwischenwert_viphj    = v(i+1,j)   + v(i,j)   ;
        zwischenwert_viph1jm1 = v(i+1,j-1) + v(i,j-1) ;
        zwischenwert_uijph    = u(i,j+1)   + u(i,j)   ;
        zwischenwert_uijmh    = u(i,j)     + u(i,j-1) ;
        zwischenwert_vimhj    = v(i,j)     + v(i-1,j) ;
        zwischenwert_uim1jph  = u(i-1,j+1) + u(i-1,j) ;

        zwischenwert_uiphj = u(i+1,j) + u(i,j)   * u(i+1,j) + u(i,j)   ;
        zwischenwert_uimhj = u(i,j)   + u(i-1,j) * u(i,j)   + u(i-1,j) ;
        zwischenwert_vijph = u(i,j+1) + u(i,j)   * u(i,j+1) + u(i,j)   ;
        zwischenwert_vijmh = u(i,j)   + u(i,j-1) * u(i,j)   + u(i,j-1) ;

        u_xx=delta_x_quadrat_divisor * (u(i+1,j) - 2 * u(i,j) + u(i-1,j));
        u_yy=delta_y_quadrat_divisor * (u(i,j+1) - 2 * u(i,j) + u(i,j-1));

        v_xx = delta_x_quadrat_divisor*(v(i+1,j) - 2 * v(i,j) + v(i-1,j));
        v_yy = delta_y_quadrat_divisor*(v(i,j+1) - 2 * v(i,j) + v(i,j-1));

        u_quadrat_x = delta_x_divisor / 4 * (zwischenwert_uiphj - zwischenwert_uimhj);
        v_quadrat_y = delta_y_divisor / 4 * (zwischenwert_vijph - zwischenwert_vijmh);

        uv_x = 1/4 * delta_x_divisor * (zwischenwert_viphj * zwischenwert_uijph - zwischenwert_vimhj   * zwischenwert_uim1jph);
        uv_y = 1/4 * delta_y_divisor * (zwischenwert_viphj * zwischenwert_uijph - zwischenwert_viphjm1 * zwischenwert_uijmh);

        u_x = delta_x_divisor * (u(i,j)-u(i-1,j));
        v_y = delta_y_divisor * (v(i,j)-v(i,j-1));

        //A(i,j)= 1/re * (u_xx(i,j)+u_yy(i,j) - uv_y(i,j) - u_quadrat_x(i,j) + gX;
        //B(i,j)= 1/re * (v_xx(i,j)+v_yy(i,j) - uv_x(i,j) - v_quadrat_y(i,j) + gY;
        
        F(i,j) = u(i,j)+ delta_t * (1/re * (u_xx(i,j)+u_yy(i,j) - uv_y(i,j) - u_quadrat_x(i,j) + gX);
        G(i,j) = v(i,j)+ delta_t * (1/re * (v_xx(i,j)+v_yy(i,j) - uv_x(i,j) - v_quadrat_y(i,j) + gY);
    }
}
