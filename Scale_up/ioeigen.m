A=dlmread('E:\�����ϸ�������\��ҵ���\MSRA\Modal Vibration\Scale up_Sparse\Scale up_Sparse\stiff.txt');
B=dlmread('E:\�����ϸ�������\��ҵ���\MSRA\Modal Vibration\Scale up_Sparse\Scale up_Sparse\mass.txt');
tic;[V,D]=eigs(A,B);toc;
X=diag(D);
fileid=fopen('E:\�����ϸ�������\��ҵ���\MSRA\Modal Vibration\Scale up_Sparse\Scale up_Sparse\eigenvalues.txt','w');
fprintf(fileid,'%d\r\n',X);
fclose(fileid);
fileid2=fopen('E:\�����ϸ�������\��ҵ���\MSRA\Modal Vibration\Scale up_Sparse\Scale up_Sparse\U_matrix.txt','w');
for i=1:size(V,2)
    p=V(:,i);
    fprintf(fileid2,'%d\r\n',p);
end
fclose(fileid2);