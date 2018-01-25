filename = 'c:\Users\amirm\Documents\Visual Studio 2017\Projects\SBS_Tools\SBS_Solver\data.json';
text = fileread(filename);
data = jsondecode(text);

Es=data.Es.real+data.Es.imag*1j;

rho=data.Rho.real+data.Rho.imag*1j;

Ep=data.Ep.real+data.Ep.imag*1j;

ax(1)=subplot(3,1,1);
plot(abs(Es(end,:)).^2-abs(Es(end,1)).^2)

ax(2)=subplot(3,1,2);
plot(abs(Ep(end,:)).^2-abs(Ep(end,1)).^2)

ax(3)=subplot(3,1,3);
plot(abs(rho(end,:)).^2-abs(rho(end,1)).^2)

