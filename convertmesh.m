function [vertices,faces] = convertmesh(D,filename,plotornot)

fid=fopen(strcat(D,'\',filename));
vertices=[];
faces=[];
for i=1:1:8
 tlines=fgetl(fid);
end
 tlines=fgetl(fid);
tlines=split(tlines);
temp=tlines{3};
NumVertex=str2num(string(temp));
tlines=fgetl(fid);
tlines=split(tlines);
temp=tlines{3};
NumTriangle=str2num(string(temp));
 while ~isequal(tlines,'[VerticesSection]')

tlines=fgetl(fid);

 end
tlines=fgetl(fid);
tlines=fgetl(fid);

for i=1:1:NumVertex
tlines=fgetl(fid);
tlines=split(tlines);
x=str2num(string(tlines{4}));
y=str2num(string(tlines{5}));
z=str2num(string(tlines{6}));
vertices=[vertices;[x y z] ];
end
for i=1:1:4
 tlines=fgetl(fid);
end
for i=1:1:NumTriangle
tlines=fgetl(fid);
tlines=split(tlines);
x=str2num(string(tlines{4}))+1;
y=str2num(string(tlines{5}))+1;
z=str2num(string(tlines{6}))+1;
faces=[faces;[x y z] ];
end

if plotornot
V_lv_epi= vertices;
F_lv_epi = faces;
c1 = [150 150 150]/255;
c1 = repmat(c1,length(V_lv_epi),1);
figure;
patch('Faces', F_lv_epi, 'Vertices', V_lv_epi, 'FaceVertexCData',c1,'FaceColor','interp');
title(filename);
ylabel('Y');
xlabel('X');
zlabel('Z')


end




fclose(fid);

end

