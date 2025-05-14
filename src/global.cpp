#include "..\include\global.hpp"

Matrix eopdata;
Matrix Cnm;
Matrix Snm;

void eop19620101(int c){
	eopdata=zeros(13,c);
	
	FILE *fp = fopen("../data/eop19620101.txt","r");
	if(fp == NULL){
		printf("Fail open eop19620101.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for(int j=1;j<=c;j++){
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&(eopdata(1,j)),
					&(eopdata(2,j)),&(eopdata(3,j)),&(eopdata(4,j)),&(eopdata(5,j)),
					&(eopdata(6,j)),&(eopdata(7,j)),&(eopdata(8,j)),&(eopdata(9,j)),
					&(eopdata(10,j)),&(eopdata(11,j)),&(eopdata(12,j)),&(eopdata(13,j)));
	}
	
	fclose(fp);
}

void GGM03S(int n){
	Cnm=zeros(n,n);
	Snm=zeros(n,n);
	
	FILE *fp = fopen("../data/GGM03S.txt","r");
	if(fp == NULL){
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	double aux;
	for(int i=1;i<=n;i++){
		for (int j=1;j<=i;j++){
			fscanf(fp,"%lf%lf%lf%lf%lf%lf",&aux,&aux,&(Cnm(i,j)),&(Snm(i,j)),&aux,&aux);
		}
	}
	
	fclose(fp);
}

Matrix PC;

void DE430Coeff(int row, int col) {
    PC = zeros(row, col);

    FILE *fp = fopen("../data/DE430Coeff.txt", "r");
    if (fp == NULL) {
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    int total_columns = 1020;

    for (int i = 1; i <= row; i++) {
        for (int j = 1; j <= col; j++) {
            fscanf(fp, "%lf", &PC(i, j));
        }
        for (int k = col + 1; k <= total_columns; k++) {
            fscanf(fp, "%lf", &aux);
        }
    }

    fclose(fp);
}
