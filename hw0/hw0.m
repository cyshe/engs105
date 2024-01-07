% Read and Plot diagonal of B
file = fopen('B_diag.dat', 'r');
formatSpec = '%f';
b_diag = fscanf(file, formatSpec);
fclose(file);

f1 = figure;
plot(b_diag, 'o');
title("Diagonal of B")
xlabel("Row")
ylabel("Value")

f2 = figure;
plot(-sqrt(b_diag), 'o');
title("Negative Square Root of Diagonal of B")
xlabel("Row")
ylabel("Value")


%%% plot third row and second subdiagonal of B
file = fopen('B.dat', 'r');
formatSpec = '%f';
B = fscanf(file, formatSpec);
fclose(file);
B = reshape(B, [6, 6]);

f3 = figure;
plot(B(3, :), 'o')
title("Third Row of B")
xlabel("Column")
ylabel("Value")

f4 = figure;
plot(diag(B, -2), 'o')
title("Second Subdiagonal of B")
xlabel("Column")
ylabel("Value")

%%%% Plot Diagonal of C
file = fopen('C_diag.dat', 'r');
formatSpec = '%f';
c_diag = fscanf(file, formatSpec);
fclose(file);

f5 = figure;
plot(c_diag, 'o');
title("Diagonal of C")
xlabel("Row")
ylabel("Value")

f6 = figure;
plot(-sqrt(c_diag), 'o');
title("Negative Square Root of Diagonal of C")
xlabel("Row")
ylabel("Value")

%%% plot third row of second subdiagonal of C
file = fopen('C.dat', 'r');
formatSpec = '%f';
C = fscanf(file, formatSpec);
fclose(file);


f7 = figure;
plot(C(1:1198), 'o')
title("Second Subdiagonal of C")
xlabel("Column")
ylabel("Value")

f8 = figure;
plot(C(1199:end), 'o')
title("Third Row of C")
xlabel("Column")
ylabel("Value")

saveas(f1,'f1.png')
saveas(f2,'f2.png')
saveas(f3,'f3.png')
saveas(f4,'f4.png')
saveas(f5,'f5.png')
saveas(f6,'f6.png')
saveas(f7,'f7.png')
saveas(f8,'f8.png')
