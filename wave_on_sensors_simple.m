function sensor_waves = wave_on_sensors_simple(cortex,  PARAMS)

    start = PARAMS.seed_vertices;
    vertices = cortex.Vertices;
    faces = cortex.Faces;
    speed = PARAMS.wave.speeds;
    duration = PARAMS.wave.duration;
    
    VertConn = tess_vertconn(vertices, faces);
    
    indn1 = find(VertConn(start,:)); 
    max_step = 9;
    num_dir = length(indn1); 
    indn = indn1;

    for n = 1:num_dir
        % find the direction vector
        IND(n,1) = start;
        ind0 = indn(n);
        IND(n,2) = ind0;
        norm0 = mean(cortex.VertNormals(indn,:));
        norm0 = norm0/norm(norm0);
        Pnorm0 = eye(3)-norm0'*norm0;
        dr0 = (cortex.Vertices(ind0,:) - cortex.Vertices(start,:));
        dr0 = dr0*Pnorm0';
        dr0 = dr0/norm(dr0);
        d = 3;
        while(d <= max_step)
            indn1 = find(VertConn(ind0,:));
            clear cs;
            for n1 = 1:length(indn1)
                dr1 = (cortex.Vertices(indn1(n1),:) - cortex.Vertices(ind0,:));
                dr1 = dr1*Pnorm0';
                dr1 = dr1/norm(dr1);

                cs(n1) = dr1*dr0';
            end;
            [csmaxval,csmaxind] = max(cs);
            ind0 = indn1(csmaxind);
            IND(n,d) = ind0;
            d = d+1;
        end
    end
        
    DIST = zeros(num_dir, (max_step-1));
    for d = 1:num_dir
        for i = 2:max_step
            DIST(d,(i-1)) = norm(vertices(IND(d,i),:)-vertices(IND(d,(i-1)),:));
        end
    end
        
    G = from_bst_get_gain_matrix(PARAMS.forward.name, PARAMS); 
    
    SR = PARAMS.sampling_rate;
    ntpoints = (SR*duration+1);
    PATH = cell(num_dir, ntpoints, length(speed));
    FM = cell(num_dir, ntpoints, length(speed));
    tstep = (1/SR);
    
    
    for s = 1:length(speed);
        l = speed(s)*tstep;
        for d = 1:num_dir
            PATH{d,1,s} = vertices(start,:);
            FM{d,1,s} = G(:,start);
            res = 0;
            v1 = 1;
            v2 = 2;
            for t = 2:ntpoints
                if l < res
                   alpha = 1-l/res;
                   PATH{d,t,s} = alpha*PATH{d,(t-1),s}+(1-alpha)*vertices(IND(d,v2),:);
                   FM{d,t,s} = alpha*FM{d,(t-1),s}+ (1-alpha)*G(:,IND(d,v2));
                   res = res-l;
                elseif l > res
                    if res == 0
                        if l < DIST(d,(v2-1))
                            alpha = 1-l/DIST(d,(v2-1));
                            PATH{d,t,s} = alpha*vertices(IND(d,v1),:)+(1-alpha)*vertices(IND(d,v2),:);
                            FM{d,t,s} = alpha*G(:,IND(d,v1)) + (1-alpha)*G(:,IND(d,v2));
                            res = DIST(d,(v2-1))-l;
                        elseif l == DIST(d,(v2-1))
                            PATH{d,t,s} = vertices(IND(d, v2),:);
                            FM{d,t,s} = G(:, IND(d, v2));
                            v1 = v1 + 1;
                            v2 = v2 + 1;
                        else
                            l2 = l-DIST(d,(v2-1));
                            v1 = v1 + 1;
                            v2 = v2 + 1;
                            alpha = 1-l2/DIST(d,(v2-1));
                            PATH{d,t,s} = alpha*vertices(IND(d,v1),:)+(1-alpha)*vertices(IND(d,v2),:);
                            FM{d,t,s} = alpha*G(:, IND(d,v1)) + (1-alpha)*G(:,IND(d,v2));
                            res = DIST(d,(v2-1))-l2;
                        end
                    else
                        l2 = l-res;
                        v1 = v1 + 1;
                        v2 = v2 + 1;
                        alpha = 1-l2/DIST(d,(v2-1));
                        PATH{d,t,s} = alpha*vertices(IND(d,v1),:)+(1-alpha)*vertices(IND(d,v2),:);
                        FM{d,t,s} = alpha*G(:,IND(d,v1))+ (1-alpha)*G(:,IND(d,v2));
                        res = DIST(d,(v2-1))-l2;
                    end
                else l == res
                    PATH{d,t,s} = vertices(IND(d, v2),:);
                    FM{d,t,s} = G(:,IND(d, v2));
                    v1 = v1 + 1;
                    v2 = v2 + 1;
                end
            end
        end
    end
            
        
%     figure
%     h = trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3));
%     set(h,'FaceAlpha',0.5);
%     s = 9;
%     hold on
%     for i = 1:num_dir
%         hplot = plot3(vertices(IND(i,:),1), vertices(IND(i,:),2), vertices(IND(i,:),3),'ko','MarkerFaceColor','r')
%     end
%     for j = 1:num_dir
%         for i = 1:ntpoints
%             hplot = plot3(PATH{j,i,s}(1), PATH{j,i,s}(2), PATH{j,i,s}(3),'ko','MarkerFaceColor','g')
%         end
%     end
    
    t = 0:(ntpoints-1);
    n = 1:ntpoints;
    
    for i = t
       wave((i+1),:) = (1 + cos(2*pi * (n - i) / ntpoints));
    end
    
%     figure
%     for i = 1:6
%     plot(t, wave(i,:))
%     hold on
%     end
       
    % waves on sensors  
    for s = 1:length(speed)
        for i = 1:num_dir
            for t = 1:ntpoints
                for k = 1:ntpoints
                    FM_s(:,k) = FM{i,k,s};
                end
                sensor_waves{i,s}(:,t) = FM_s*wave(t,:)';
            end
        end
    end
    
    % spherical wave
    for s = 1:length(speed)
        sensor_waves{(num_dir+1),s} = zeros(size(sensor_waves{1,1}));
        for t = 1:ntpoints
            for i = 1:num_dir
                sensor_waves{(num_dir+1),s}(:,t) = sensor_waves{(num_dir+1),s}(:,t)+sensor_waves{i,s}(:,t);
            end
        end
    end
    
    
    
end





