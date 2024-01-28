
function u = FGLSE(Img,guidance,rt,Initial_LSF,Iternum,sigma,lamda1_loc,lamda2_loc,Kernel_diffusion,epsilon)
    Img = double(Img(:,:,1));
    [row,col] = size(Img);
    
    delt_t = 0.1; 
    lamda1_gol = 1.05; lamda2_gol = 1.00;
    Kernel = fspecial('gaussian',round(2*sigma)*2+1,sigma);
    u = Initial_LSF;

    for i = 1:Iternum   
        Drc_u = (epsilon/pi)./(epsilon^2.+u.^2); 

        % local temp
        [f1,f2,local] = LRCV(Img,u,epsilon,Kernel,lamda1_loc,lamda2_loc);

        % global temp    
        dx = 1./(0.1 + abs(f1-f2));  
        golbal_fit = CV_fit(guidance,u,epsilon,lamda1_gol,lamda2_gol);
        golbal = rt(i).* golbal_fit.*dx;

        % Iterating
        u0 = u;
        u = u + delt_t* Drc_u.*(golbal +local);
        u = imfilter(u,Kernel_diffusion,'replicate','same');
 
        relative_error = (double(u0>0)-double(u>0))/(sum(sum(double(u0>0)))+0.001);;
        if relative_error <0.00001  
            break;
        end    
    end
end
%% º¯ÊýºÏ¼¯

function golbal = CV_fit(Img,u,epsilon,lamda1,lamda2)
    H_u = 0.5*(1+(2/pi)*atan(u/epsilon));
    p1 = sum(sum(Img.*H_u));
    p2 = sum(sum(H_u));
    p3 = sum(sum(Img));
    [m,n] = size(Img);
    p4 = m*n;
    c1 = p1/p2;
    c2 = (p3-p1) / (p4-p2);
    golbal = -lamda1*(Img - c1).^2 + lamda2*(Img-c2).^2;
end


 function [f1,f2,dataTerm] = LRCV(Img,u,epsilon,Kernel,lamda_1,lamda_2)
    H = 0.5*(1+(2/pi)*atan(u./epsilon));
    [f1,f2] = Local_Avr(Img,H,Kernel);
    dataTerm = -lamda_1*(Img - f1).^2 + lamda_2*(Img - f2).^2;

    function [f1,f2] = Local_Avr(I,H,K)
        f1 = conv2(I.*H,K,'same');
        c1 = conv2(H,K,'same');
        f1 = f1./c1;
        f2 = conv2(I.*(1-H),K,'same');
        c2 = conv2(1-H,K,'same');
        f2 = f2./c2;
    end
 end
 



            
