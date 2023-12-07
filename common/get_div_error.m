function gauss_err = get_div_error(Ex,Ey,rhoOverSigma1,dx,dy)
    E_x = diff(Ex,1,1)/dx;
    E_y = diff(Ey,1,2)/dy;
    divE = E_x + E_y;

    gauss_err = rhoOverSigma1(2:end-1,2:end-1) - divE; % assuming rho has already been divided by sigma_1
end