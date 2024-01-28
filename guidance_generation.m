
function guidance = guidance_generation(Img,tao,n)
    Img = double(Img(:,:,1));
    [dx,dy] = gradient(Img);
    dImg_org = max(abs(dx),abs(dy)); 
    M = max(max(dImg_org));
    h =  (1-exp(-((dImg_org/M))*tao));
    fmin = noise_if(Img);
    v = h.*fmin;

    Img_Dv_mag = GL(Img,n,v);
    guidance = Img - Img_Dv_mag;
    guidance = double(guidance>=0).*guidance;
   
end


function diff = noise_if(Img)
    [row,col] = size(Img);
    diff = ones(row,col);
    Img_max = ordfilt2(Img,9,ones(3,3));
    Img_min = ordfilt2(Img,1,ones(3,3));
    diff_max = Img_max-Img_min;
    for i = 2:(row-1)
        for j = 2:(col-1)
            if diff_max(i,j) == 0
                diff(i,j) = 1;
             else
                temp = [abs(Img(i+1,j)-Img(i-1,j)),abs(Img(i+1,j+1)-Img(i-1,j-1)),abs(Img(i,j+1)-Img(i,j-1)),abs(Img(i+1,j-1)-Img(i-1,j+1))];
                diff(i,j) = max(temp)/(diff_max(i,j));
            end
        end
    end
    diff = ordfilt2(diff,1,ones(3,3),'symmetric');
end

function Mag_Dv = GL(Img,m,v)
% m是所用到前m个像素值，v是自适应分数阶矩阵
    [row,col] = size(Img);
    % 利用m对Img扩充
    Img_ext = [Img(:,m+1:-1:2),Img,Img(:,col-1:-1:col-m)];
    Img_ext = [Img_ext(m+1:-1:2,:);Img_ext;Img_ext(row-1:-1:row-m,:)];
    Dv_x_1 = Img; Dv_y_1 = Img; Dv_x_2 = Img; Dv_y_2 = Img;
    for i  = 1:row
        for j = 1:col
            % 形成系数
            w(1) = 1;
            for k = 2:m+1
                w(k) = w(k-1)*(-v(i,j)+(k-2))/(k-1);
                Dv_x_1(i,j) = Dv_x_1(i,j) + w(k)*Img_ext(i+m,j+m-k+1);
                Dv_x_2(i,j) = Dv_x_2(i,j) + w(k)*Img_ext(i+m,j+m+k-1);
                
                Dv_y_1(i,j) = Dv_y_1(i,j) + w(k)*Img_ext(i+m-k+1,j+m);  
                Dv_y_2(i,j) = Dv_y_2(i,j) + w(k)*Img_ext(i+m+k-1,j+m);  
            end
        end
    end
    Dv_x = 0.5*(Dv_x_1+Dv_x_2);    Dv_y = 0.5*(Dv_y_1+Dv_y_2);
    Mag_Dv = (Dv_x.^2+Dv_y.^2).^(1/2); 
end
